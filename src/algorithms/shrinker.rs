//! Shrink the given graph
use crate::graphs::Graph;
use crate::graphs::pt_graph::{CarryType, EdgeIndex, EdgeType, EmptyEdge, NodeIndex, PtGraph};
use crate::logw;

use petgraph::Direction::{Incoming, Outgoing};
use petgraph::visit::EdgeRef;

// use std::process::exit;

/// Mark graph as shrinkable.
pub trait Shrinkable {
    /// Edge index associated with collection.
    type EdgeIdx;
    /// Node index associated with collection.
    type NodeIdx;
    /// Shrink graph.
    ///
    /// This operation should shrink all straight paths
    /// It is assumed that after shrinking graph will not have any nodes
    /// connected in this way: s -> x -> ... -> t
    fn shrink(&mut self) -> bool;

    /// Shrink one single path. This method assumes that `base_edge` argument points
    /// to a valid edge, which target has a single outgoing edge.
    /// Returns index of the shrinked path represented by edge.
    fn shrink_single_path(&mut self, start_node: Self::NodeIdx, mid_node : Self::NodeIdx, ambnodes : &<PtGraph as Graph>::AmbiguousNodes, currtype : EdgeType);

    /// This helper function modifies the edges accordingly to have shrunken nodes with
    /// internal edges.
    fn modify_edges_when_shrinking(&mut self, base_node: NodeIndex, prev_node : NodeIndex, internal_edge_ty : EdgeType, in_edge_ind : EdgeIndex, in_edge_ty : EdgeType);

}

impl Shrinkable for PtGraph {
    type EdgeIdx = EdgeIndex;
    type NodeIdx = NodeIndex;

    fn shrink(&mut self) -> bool {

        // Shrinkage here means to only find consecutive nodes, w/o bifurcations

        let mut dididoanything = false;
        loop {
            let mut dididoanythingnow = false;
            let ambnodes = self.get_ambiguous_nodes_bi(); // Just in case we hadn't got them yet
            logw(format!("Starting shrinking the graph with {} nodes and {} edges, beginning from {} ambiguous nodes",
                  self.graph.node_count(),
                  self.graph.edge_count(),
                  ambnodes.len()).as_str(), Some("trace"));
            for an in ambnodes.iter() {
                // println!("it");
                // We need to see the forward edges from here

                if !self.graph.contains_node(*an) {
                    continue;
                }

                let neigh = self.get_good_neighbours_bi(*an);
                let conns = neigh.len();

                if conns == 1 {
                    // println!("1 conn");
                    // This means that this is the outermost k-mer of a dead-end, that might also have self-loops.
                    // I.e. this k-mer itself, if no self-loops are present, can be susceptible of being shrunk.

                    // Let's check for self-loops first:
                    if self.node_has_self_loops(*an) {
                        continue;
                    }

                    let tmpty = neigh[0].1.get_from_and_to().1;
                    let outn = self.out_neighbours_bi(neigh[0].0, tmpty);
                    if outn.len() == 1 && (outn[0].0 == *an || ambnodes.contains(&outn[0].0)) {
                        continue;
                    } else if outn.len() <= 1 && self.in_neighbours_bi(neigh[0].0, tmpty).len() == 1 {
                        // println!("yes1");
                        self.shrink_single_path(*an, neigh[0].0, &ambnodes, neigh[0].1);
                        dididoanything = true;
                        dididoanythingnow = true;
                    }

                } else {
                    // println!("various conn");
                    // This situation means that apart from possible self-loops, we have various valid connections. This
                    // k-mer is not valid itself for shrinkage, but one of those with which is connected might be.

                    // println!("> number of true neighbours: {:?}", neigh.len(),);
                    for n in neigh {
                        // println!("\t- neighid {:?}, type of edge that connects an with it: {:?}", n.0, n.1);

                        if ambnodes.contains(&n.0) || n.0 == *an {
                            // println!("\t- It's an ambiguous node!");
                            continue;
                        } else {
                            // println!("\t- It's NOT an ambiguous node!");
                            let tmpty = n.1.get_from_and_to().1;
                            let outn = self.out_neighbours_bi(n.0, tmpty);

                            if outn.len() == 0 {
                                continue;
                            }

                            if outn.len() == 1 && outn[0].0 != n.0 && outn[0].0 != *an &&
                                (!ambnodes.contains(&outn[0].0) || self.get_good_neighbours_bi(*an).len() == 1) { // This last condition allow us to catch
                                                                                            // cases where from a junction you try to shrink a path that ends.
                                let incn = self.in_neighbours_bi(n.0, tmpty);
                                if incn.len() == 1 {
                                    // println!("yes2");
                                    self.shrink_single_path(n.0, outn[0].0, &ambnodes, outn[0].1);
                                    dididoanything = true;
                                    dididoanythingnow = true;
                                }
                            }
                        }
                    }
                }
            }
            if !dididoanythingnow {
                break;
            }
        }


        log::info!("Shrinking ended. Shrunk graph has {} nodes and {} edges",
              self.graph.node_count(),
              self.graph.edge_count());

        dididoanything
    }



    #[inline]
    fn modify_edges_when_shrinking(&mut self, base_node: NodeIndex, prev_node : NodeIndex, internal_edge_ty : EdgeType, in_edge_ind : EdgeIndex, in_edge_ty : EdgeType) {

        log::trace!("base_node: {:?} prev_node: {:?} internal_edge_ty: {:?} in_edge_ind: {:?} in_edge_ty: {:?}",
            base_node, prev_node, internal_edge_ty, in_edge_ind, in_edge_ty,
        );
        log::trace!("edges coming to the prev_node from the base_node: {:?}", self.graph.edges_connecting(base_node, prev_node).map(|x| x.id()).collect::<Vec<_>>());
        log::trace!("edges coming to the base_node from the prev_node: {:?}", self.graph.edges_connecting(prev_node, base_node).map(|x| x.id()).collect::<Vec<_>>());

        match internal_edge_ty {
            EdgeType::MinToMax => {
                match in_edge_ty {
                    EdgeType::MinToMin  => {
                        self.graph.edge_weight_mut(in_edge_ind).unwrap().t = EdgeType::MinToMax;
                        self.graph.edge_weight_mut(self.graph.edges_connecting(base_node, prev_node).next().unwrap().id()).unwrap().t = EdgeType::MinToMax;
                    },
                    EdgeType::MaxToMin  => {
                        self.graph.edge_weight_mut(in_edge_ind).unwrap().t = EdgeType::MaxToMax;
                        self.graph.edge_weight_mut(self.graph.edges_connecting(base_node, prev_node).next().unwrap().id()).unwrap().t = EdgeType::MinToMin;
                    },
                    _                   => panic!("Not expected edge type"),
                }
                self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MaxToMax);
            },
            EdgeType::MaxToMin => {
                match in_edge_ty {
                    EdgeType::MaxToMax  => {
                        self.graph.edge_weight_mut(in_edge_ind).unwrap().t = EdgeType::MaxToMin;
                        self.graph.edge_weight_mut(self.graph.edges_connecting(base_node, prev_node).next().unwrap().id()).unwrap().t = EdgeType::MaxToMin;
                    },
                    EdgeType::MinToMax  => {
                        self.graph.edge_weight_mut(in_edge_ind).unwrap().t = EdgeType::MinToMin;
                        self.graph.edge_weight_mut(self.graph.edges_connecting(base_node, prev_node).next().unwrap().id()).unwrap().t = EdgeType::MaxToMax;
                    },
                    _                   => panic!("Not expected edge type"),
                }
                self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MinToMin);
            },
            _ => panic!("Value not expected"),
        }
    }



    #[inline]
    fn shrink_single_path(&mut self, base_node: NodeIndex, mut next_node: NodeIndex,
                          ambnodes : &<PtGraph as Graph>::AmbiguousNodes, mut curredge : EdgeType) {

        let mut countsformean : Vec<u16> = vec![self.graph.node_weight(base_node).unwrap().counts];

        let (initty, mut currtype) = curredge.get_from_and_to();

        log::trace!("Starting shrinkage!");
        if self.out_degree_bi(next_node, currtype) == 0 {
            // Only the base_node + next_node shrinkage is possible
            log::trace!("No one to continue already from the next_node, before looping");

            countsformean.push(self.graph.node_weight(next_node).unwrap().counts);
            let next_base_weight = self.graph.remove_node(next_node).unwrap();
            let ind = self.in_neighbours_bi(base_node, initty);

            // self.graph.node_weight_mut(base_node).unwrap().merge(&next_base_weight);             // OLD
            self.graph.node_weight_mut(base_node).unwrap().merge(&next_base_weight, curredge);      // NEW
            self.graph.node_weight_mut(base_node).unwrap().set_mean_counts(&countsformean);


            log::trace!("base_node {:?} next_node {:?} curredge {:?} initty {:?} currtype {:?} ind {:?}",
                        base_node, next_node, curredge, initty, currtype, ind);


            // ======================================= CHANGING INTERNAL EDGES IF NEEDED BEGIN

            self.graph.node_weight_mut(base_node).unwrap().invert_if_needed(curredge);

            if curredge.is_direct() {
                self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(curredge);

            } else if ind.len() == 0 {
                // NOTE: now, this is set to MintoMin by default. IT IS A LIE
                self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MinToMin);

            } else {
                // The base_node is connected to other nodes, previously
                // AND the edge is not direct, it's crossed. Therefore, we
                // need to change those edges connecting to the previous
                // node.

                // NOTE: this CHANGES the edges between nodes, and makes that
                // them are FALSE

                let theedges : Vec<_> = self.graph.edges_connecting(ind[0].0, base_node).map(|e| e.id()).collect();
                if theedges.len() > 1 {
                    panic!("More than one linking outgoing edge, this should not happen unless there are multiple connections to the same node.");
                }

                log::trace!("Modifying edges with a non-direct edge at the beginning.");
                self.modify_edges_when_shrinking(base_node, ind[0].0, curredge, theedges[0], ind[0].1);

            }
            // ======================================= CHANGING INTERNAL EDGES IF NEEDED END


            // =================== DEBUG
            if self.in_degree(base_node) != 0 || self.out_degree(base_node) != 0 {
                panic!("Nope! 1");
            }
            // =================== DEBUG

            return;
        }

        // As next_node has exactly one outgoing neighbour, we can start the main loop
        // =================== DEBUG
        if self.in_degree(next_node) != 2 || self.out_degree(next_node) != 2 {
            panic!("Nope!");
        }
        // =================== DEBUG

        // println!("Looping!");
        loop {
            // println!("\t- New iteration");
            // The out_degree MUST be 1 when you arrive here, and also the in_degree MUST be 1 (of the next_node)
            // Note that these conditions here force that:
            // THIS FOLLOWING COMMENT LINES ARE OLD, TODO REVIEW
            //     - A perfect loop of hashes will get shrunk to two nodes connected by two directed edges IF you
            // feed such to this function, a thing that at the time of writing this, is not implemented.
            //     - Self-loops might still exist, but will only be possible if a kmer has, as forward neighbour,
            // itself (always modulo collisions). Thus, these are much more unusual.
            //     - (TODO) situations with two nodes connected by only one directed edge, in which the origin node
            // of such edge has zero, one, or more than one incoming edges and only one outgoing edge, and in which
            // the destination node has only one incoming edge and zero, one, or more outgoing edges, will NOT get
            // shrunk under the current code.

            // Let's get our potential next node!
            let prospective_node = self.out_neighbours_bi(next_node, currtype)[0];

            // First, avoid self-loops. This could be changed in the future.
            if prospective_node.0 == base_node {
                panic!("This should not happen");
            }

            countsformean.push(self.graph.node_weight(next_node).unwrap().counts);
            let next_base_weight = self.graph.remove_node(next_node).unwrap();
            let prospfromandto = prospective_node.1.get_from_and_to();
            currtype = prospfromandto.0;

            self.graph.node_weight_mut(base_node).unwrap().merge(&next_base_weight, curredge);


            if self.out_degree_bi(prospective_node.0, prospfromandto.1) == 0 && self.in_degree_bi(prospective_node.0, prospfromandto.1) == 0 {
                // Prospective_node is an ambiguous node, but only because it is an external where we can finish. Let's do that.
                // println!("Prospective is ambiguous and does not have any other continuation: we'll finish.");
                countsformean.push(self.graph.node_weight(prospective_node.0).unwrap().counts);
                let next_base_weight = self.graph.remove_node(prospective_node.0).unwrap();
                let nw = self.graph.node_weight_mut(base_node).unwrap();

                curredge = EdgeType::from_carrytypes(initty, prospfromandto.1); // Because we include the prosp node
                nw.merge(&next_base_weight, curredge);

                nw.set_mean_counts(&countsformean);

                // ======================================= CHANGING INTERNAL EDGES IF NEEDED BEGIN
                let ind = self.in_neighbours_bi(base_node, initty);
                // println!("Length of incoming neighbours {:?}", ind.len());

                self.graph.node_weight_mut(base_node).unwrap().invert_if_needed(curredge);

                if curredge.is_direct(){
                    self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(curredge);
                } else if ind.len() == 0 {
                    // NOTE: now, this is set to MintoMin by default. IT IS A LIE
                    self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MinToMin);
                } else {
                    // The base_node is connected to other nodes, previously
                    // AND the edge is not direct, it's crossed. Therefore, we
                    // need to change those edges connecting to the previous
                    // node.

                    // NOTE: this CHANGES the edges between nodes, and make than
                    // ones of them are FALSE

                    let theedges : Vec<_> = self.graph.edges_connecting(ind[0].0, base_node).map(|e| e.id()).collect();
                    if theedges.len() > 1 {
                        panic!("More than one linking outgoing edge, this should not happen unless there are multiple connections to the same node.");
                    }

                    log::trace!("Modifying edges with a non-direct edge in the loop to an ambiguous node that is an external");
                    self.modify_edges_when_shrinking(base_node, ind[0].0, curredge, theedges[0], ind[0].1);
                }
                // ======================================= CHANGING INTERNAL EDGES IF NEEDED END

                // =================== DEBUG
                if self.out_degree(base_node) != 1 && ind.len() != 0 {
                    println!("{:?} {:?}", self.out_degree(base_node), self.in_degree(base_node));
                    panic!("EY2!");
                }
                // =================== DEBUG

                return;

            } else if ambnodes.contains(&prospective_node.0) {
                // We cannot add prospective_node, because it is an ambiguous node. So as we have already added the
                // info from next_node to base_node, we just have to set up the internal edge and link properly
                // base_node to the prospective_node.

                self.graph.node_weight_mut(base_node).unwrap().set_mean_counts(&countsformean);


                // ======================================= CHANGING INTERNAL EDGES IF NEEDED BEGIN
                let ind = self.in_neighbours_bi(base_node, initty);
                curredge = EdgeType::from_carrytypes(initty, currtype);

                let newoutedge : EdgeType;

                self.graph.node_weight_mut(base_node).unwrap().invert_if_needed(curredge);

                if curredge.is_direct() {
                    self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(curredge);
                    newoutedge = prospective_node.1;
                    // println!("DIRECT");
                } else if ind.len() == 0 {
                    // We need to change the FOLLOWING edges, not the previous ones!
                    // println!("indlenzero");
                    match curredge {
                        EdgeType::MinToMax => {
                            newoutedge = EdgeType::from_carrytypes(CarryType::Min, prospective_node.1.get_from_and_to().1);
                            self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MinToMin);
                        },
                        EdgeType::MaxToMin => {
                            newoutedge = EdgeType::from_carrytypes(CarryType::Max, prospective_node.1.get_from_and_to().1);
                            self.graph.node_weight_mut(base_node).unwrap().set_internal_edge(EdgeType::MaxToMax);
                        },
                        _ => panic!("Value not expected"),
                    }
                } else {
                    // println!("indlenNOTzero");
                    newoutedge = prospective_node.1;
                    // The base_node is connected to other nodes, previously
                    // AND the edge is not direct, it's crossed. Therefore, we
                    // need to change those edges connecting to the previous
                    // node.

                    // NOTE: this CHANGES the edges between nodes, and make that
                    // a pair of those are FALSE

                    let theedges : Vec<_> = self.graph.edges_connecting(ind[0].0, base_node).map(|e| e.id()).collect();
                    if theedges.len() > 1 {
                        panic!("More than one linking outgoing edge, this should not happen unless there are multiple connections to the same node.");
                    }

                    log::trace!("Modifying edges with a non-direct edge in the loop when reaching another ambiguous node that is not an external.");
                    self.modify_edges_when_shrinking(base_node, ind[0].0, curredge, theedges[0], ind[0].1);
                }

                self.graph.add_edge(base_node, prospective_node.0, EmptyEdge{t : newoutedge});
                self.graph.add_edge(prospective_node.0, base_node, EmptyEdge{t : newoutedge.rev()});
                // ======================================= CHANGING INTERNAL EDGES IF NEEDED END

                // =================== DEBUG
                if self.out_degree_min(base_node) != 1 && ind.len() != 0 || self.out_degree(base_node) != 1 && ind.len() == 0{
                    println!("Base node: {:?}, internal type: {:?}", base_node, self.graph.node_weight(base_node).unwrap().innerdir.unwrap());
                    println!("Incoming node to base_node from before: {:?}, prospective_node: {:?}", ind[0].0, prospective_node.0);
                    println!("OUTGOING");
                    for e in self.graph.edges_directed(base_node, Outgoing) {
                        println!("- Target: {:?} Type: {:?}", e.target(), e.weight().t);
                    }
                    println!("INCOMING");
                    for e in self.graph.edges_directed(base_node, Incoming) {
                        println!("- Source: {:?} Type: {:?}", e.source(), e.weight().t);
                    }

                    panic!("EY1!");
                }
                // =================== DEBUG

                return;

            } else {
                // If we are here, that means we can continue shrinking!
                // =================== DEBUG
                if self.in_degree(prospective_node.0) != 1 || self.out_degree(prospective_node.0) != 1 {
                    println!("{:?} {:?}", self.in_degree(prospective_node.0), self.out_degree(prospective_node.0));
                    panic!("Nope!");
                }
                // println!("Continuing! in_deg {} out_deg {}", self.in_degree(prospective_node.0), self.out_degree(prospective_node.0));
                // =================== DEBUG

                // Updating next_node and currtype
                next_node = prospective_node.0;
                currtype  = prospective_node.1.get_from_and_to().1;
                curredge  = EdgeType::from_carrytypes(initty, currtype);
            }
        }
    }
}
