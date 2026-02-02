//! Corrects parts of the provided graph, if needed
use crate::graphs::pt_graph::{
    CarryType, EdgeIndex, EdgeType, EmptyEdge, NodeIndex, NodeStruct, PtGraph,
};
use crate::graphs::Graph;
use crate::{logw, EdgeWeight};

use petgraph::visit::EdgeRef;
use petgraph::Direction::{Incoming, Outgoing};
use std::{
    cmp::{max, min},
    collections::BTreeSet,
    vec::Drain,
};

// use std::process::exit;

/// Mark graph as shrinkable.
pub trait Correctable {
    /// Edge index associated with collection.
    type EdgeIdx;

    /// Node index associated with collection.
    type NodeIdx;

    /// Remove edges with weight below threshold.
    fn remove_weak_nodes(&mut self, threshold: EdgeWeight);

    /// Remove edges that are self-loops, i.e. those whose source and destination nodes are the same.
    fn remove_self_loops(&mut self);

    /// Solve bubbles from the graph
    fn correct_bubbles(&mut self) -> bool;

    /// Remove all input and output dead paths
    fn remove_dead_paths(&mut self) -> bool;

    /// Find and remove all links that are impossible in bi-directed de Bruijn graphs derived from DNA sequences.
    fn remove_conflictive_links(&mut self) -> bool;
}

impl Correctable for PtGraph {
    type EdgeIdx = EdgeIndex;

    type NodeIdx = NodeIndex;

    fn remove_weak_nodes(&mut self, threshold: EdgeWeight) {
        self.graph
            .retain_nodes(|g, n| g.node_weight(n).unwrap().counts >= threshold);
    }

    fn remove_self_loops(&mut self) {
        self.graph.retain_edges(|g, e| {
            let (n1, n2) = g.edge_endpoints(e).unwrap();
            n1 != n2
        });
    }

    fn correct_bubbles(&mut self) -> bool {
        let mut dididoanything = false;

        logw("Starting resolution of standard bubbles", Some("info"));
        let bubbles = self
            .graph
            .node_indices()
            .filter(|n| self.out_degree(*n) == 3)
            .filter(|n| {
                // We can directly check only once the bubbles (and not do it twice) by fixing
                // the carrytype from which we want to continue
                let mut vmin = Vec::with_capacity(2);

                for e in self.graph.edges_directed(*n, Outgoing) {
                    if e.weight().t.get_from_and_to().0 == CarryType::Min {
                        if vmin.len() == 2 {
                            return false;
                        } else {
                            vmin.push(e.id());
                        }
                    }
                }

                // We have two outgoing from min/max and one outgoing from max/min.
                // We check whether this is actually a bubble or not.
                if vmin.len() == 2 {
                    check_bubble_structure(self, *n, vmin)
                } else {
                    false
                }
            })
            .collect::<BTreeSet<NodeIndex>>();

        if bubbles.is_empty() {
            return false;
        } else {
            logw(
                format!(
                    "Found {:?} potential bubbles (they might be less). Starting to collapse them ",
                    bubbles.len()
                )
                .as_str(),
                Some("trace"),
            );
            for n in bubbles {
                if self.graph.contains_node(n) {
                    let tmpb = collapse_bubble(self, n);
                    if tmpb {
                        dididoanything = true;
                    }
                }
            }
        }

        logw(
            format!(
                "Bubble correction ended. Corrected graph has {} nodes and {} edges",
                self.graph.node_count(),
                self.graph.edge_count()
            )
            .as_str(),
            Some("info"),
        );

        dididoanything
    }

    fn remove_dead_paths(&mut self) -> bool {
        logw(
            format!(
                "Before pruning: {} nodes and {} edges",
                self.graph.node_count(),
                self.graph.edge_count()
            )
            .as_str(),
            Some("info"),
        );

        let mut dididoanything = false;
        logw(format!("Starting graph pruning. Graph has {} externals, {} alone nodes, the remaining are internal.",
            self.externals_bi().len(),
            self.graph.node_indices().filter(|n| self.out_degree(*n) == 0 && self.in_degree(*n) == 0).count()).as_str(), Some("info"));

        let mut to_remove: Vec<NodeIndex> = vec![];
        loop {
            let mut path_check_vec = vec![];
            let externals: Vec<_> = self
                .externals_bi()
                .into_iter()
                .filter(|n| self.out_degree(*n) == 1)
                .collect();

            // TODO: go back to "Graph" alone, for speed

            logw(
                format!("Detected {} externals", externals.len()).as_str(),
                Some("trace"),
            );

            for v in externals {
                check_dead_path(
                    self,
                    v,
                    &mut path_check_vec,
                    self.k,
                    self.graph
                        .edges_directed(v, petgraph::EdgeDirection::Outgoing)
                        .next()
                        .unwrap()
                        .weight()
                        .t,
                );
                if !path_check_vec.is_empty() {
                    dididoanything = true;
                    to_remove.append(&mut path_check_vec);
                }
            }

            // if there are no dead paths left pruning is done
            if to_remove.is_empty() {
                logw("Graph is pruned.", Some("info"));
                return dididoanything;
            }

            // reverse sort edge indices such that removal won't cause any troubles with swapped
            // edge indices (see `petgraph`'s explanation of `remove_edge`)
            // NOTE: this might be useful for the change to Graph backend, but now we don't need it.
            // to_remove.sort_by(|a, b| b.cmp(a));
            remove_paths(self, to_remove.drain(..));
        }
    }

    fn remove_conflictive_links(&mut self) -> bool {
        // Conflicting links: k-mers connected to others with a relatively larger count
        // TODO: apply
        false

        // OLD-correct bad links (only one edge instead of two, or bad types) BEGIN
        // let mut dididoanything = false;
        // log::info!("Starting removal of conflicting links. Graph has {} sources, {} sinks, {} alone nodes, the remaining are internal.",
        //     self.graph.externals(EdgeDirection::Incoming).count(), self.graph.externals(EdgeDirection::Outgoing).count(),
        //     self.graph.node_indices().filter(|n| self.out_degree(*n) == 0 && self.in_degree(*n) == 0).count());
        //
        // log::trace!("Detected {} sources and {} sinks", self.graph.externals(EdgeDirection::Incoming).count(), self.graph.externals(EdgeDirection::Outgoing).count());
        //
        // let susceptiblenodes : Vec<NodeIndex> = self.graph.node_indices().filter(|&n| {
        //
        //     self.graph.neighbors_directed(n, EdgeDirection::Incoming).count() > 1 || self.graph.neighbors_directed(n, EdgeDirection::Outgoing).count() > 1
        //
        // }).collect::<Vec<NodeIndex>>();
        //
        // for v in susceptiblenodes {
        //     let incneigh = self.graph.neighbors_directed(v, EdgeDirection::Incoming).collect::<Vec<NodeIndex>>();
        //     let outneigh = self.graph.neighbors_directed(v, EdgeDirection::Outgoing).collect::<Vec<NodeIndex>>();
        //
        //
        //     let suscnodecs = self.graph.node_weight(v).unwrap().counts;
        //     // println!("\nSusc. node counts: {}", suscnodecs);
        //     for incn in incneigh.iter() {
        //         let coc = self.graph.node_weight(*incn).unwrap().counts as f32/suscnodecs as f32;
        //         // if coc < 0.1 {
        //         //     println!("\tInc. node count: {}\t{:.2}", self.graph.node_weight(*incn).unwrap().counts, coc);
        //         // }
        //     }
        //
        //     for outn in outneigh.iter() {
        //         let coc = self.graph.node_weight(*outn).unwrap().counts as f32/suscnodecs as f32;
        //         // if coc < 0.1 {
        //         //     println!("\tOut. node count: {}\t{:.2}", self.graph.node_weight(*outn).unwrap().counts, coc);
        //         // }
        //     }
        // }
        //
        // return dididoanything

        // OLD END

        // TEST TEMP for creating function for removing conflictive links, but it is better to do it while shrinking the graph the first time
        // let mut checkednodes : BTreeSet<NodeIndex> = BTreeSet::new();
        // let mut dididoanything = false;
        //
        // for ni in self.ptgraph.node_indices() {
        //     if check_connections_and_remove(&mut self.ptgraph, ni, self.in_neighbours_min(ni), self.out_neighbours_min(ni), &checkednodes) {
        //         dididoanything = true;
        //     }
        //     checkednodes.insert(ni);
        // }
        //
        // return dididoanything;
    }
}

// TEST TEMP for creating function for removing conflictive links, but it is better to do it while shrinking the graph the first time
// fn check_connections_and_remove(ptgraph: &mut PtGraph, ni : NodeIndex, inneighs : Vec<(NodeIndex, EdgeType)>, outneighs : Vec<(NodeIndex, EdgeType)>, nodeset : &BTreeSet<NodeIndex>) -> bool {
//     let mut didierasedanything = false;
//     let nodecounts = ptgraph.node_weight(ni).unwrap().counts as f32;
//
//     for (neigh, edge) in inneighs {
//         if ptgraph.node_weight(neigh).unwrap().counts as f32 >
//     }
//
//
//
//
//     return didierasedanything;
// }

/// Checks whether the candidate area can be a good bubble for error correction.
fn check_bubble_structure(ptgraph: &PtGraph, startn: NodeIndex, invec: Vec<EdgeIndex>) -> bool {
    let mut midnodes = Vec::with_capacity(2);
    let mut midcts = Vec::with_capacity(2);

    for e in invec {
        midnodes.push(ptgraph.graph.edge_endpoints(e).unwrap().1);
        midcts.push(ptgraph.graph.edge_weight(e).unwrap().t.get_from_and_to().1);
    }

    // We need to check that the two intermediate nodes are different
    if midnodes[0] == midnodes[1] || midnodes[0] == startn || midnodes[1] == startn {
        return false;
    }

    // Now, how many neighbours do we have from the middle nodes? It should be only one, and the same
    // and also, only one incoming neighbour!
    let tmpv0 = ptgraph.out_neighbours_bi(midnodes[0], midcts[0]);
    if tmpv0.len() != 1 || ptgraph.in_neighbours_bi(midnodes[0], midcts[0]).len() != 1 {
        return false;
    }
    let outnode = tmpv0[0].0;
    let tmpv1 = ptgraph.out_neighbours_bi(midnodes[1], midcts[1]);
    let outct = tmpv0[0].1.get_from_and_to().1;

    if (tmpv1.len() != 1)
        || (tmpv1[0].0 != outnode)
        || (outct != tmpv1[0].1.get_from_and_to().1)
        || (ptgraph.in_neighbours_bi(midnodes[1], midcts[1]).len() != 1)
    {
        return false;
    }

    // and now we check that we go to a new node that is different from outnode (i.e. no self-edges/loops)
    let tmpv3 = ptgraph.out_neighbours_bi(outnode, outct);
    if tmpv3.len() != 1
        || tmpv3[0].0 == startn
        || ptgraph.in_neighbours_bi(outnode, outct).len() != 2
    {
        return false;
    }

    true
}

/// This function collapses standard bubbles depending on the number of counts (very naive)
fn collapse_bubble(ptgraph: &mut PtGraph, startn: NodeIndex) -> bool {
    // Remember: we always start with the Min being the origin of the two outward edges
    // let prevconn = ptgraph.in_neighbours_min(startn)[0];
    let midconns = ptgraph.out_neighbours_min(startn);
    if midconns.len() != 2
        || ptgraph
            .out_neighbours_bi(midconns[0].0, midconns[0].1.get_from_and_to().1)
            .len()
            != 1
    {
        return false;
    }
    let chosennode: usize;
    let node0w = ptgraph.graph.node_weight(midconns[0].0).unwrap();
    let node1w = ptgraph.graph.node_weight(midconns[1].0).unwrap();
    let savedmidw: NodeStruct;

    let count_threshold = min(
        (0.1 * max(node0w.counts, node1w.counts) as f32).round() as u16,
        1,
    ); // Inspired by Skesa

    // println!("here1");

    if node0w.counts < count_threshold {
        chosennode = 1;
        savedmidw = node1w.clone();
        // println!("here3");
    } else if node1w.counts < count_threshold {
        chosennode = 0;
        savedmidw = node0w.clone();
        // println!("here4");
    } else {
        // println!("here2");
        // We have not managed to exclude one of the two options. What do we do now?

        if ((node0w.abs_ind.len() - node1w.abs_ind.len()) as i32).abs() as f32
            / (max(node0w.abs_ind.len(), node1w.abs_ind.len()) as f32)
            > 0.025
        {
            // There is a large difference in the number of k-mers between the options: this could be a coincidence link between
            // two locations far away in the genome. Let's compare with the counts of the start and end nodes.
            // println!("here5");

            let startn_counts = ptgraph.graph.node_weight(startn).unwrap().counts;

            // println!("preasdf");
            let endn_counts = ptgraph
                .graph
                .node_weight(
                    ptgraph.out_neighbours_bi(midconns[0].0, midconns[0].1.get_from_and_to().1)[0]
                        .0,
                )
                .unwrap()
                .counts;
            // println!("asdfadsf");
            let average_surrounding_counts =
                ((startn_counts + endn_counts) as f32 / 2.0).round() as u16;

            let rel_diff_0 = ((node0w.counts - average_surrounding_counts) as i32).abs() as f32
                / (average_surrounding_counts as f32);
            let rel_diff_1 = ((node1w.counts - average_surrounding_counts) as i32).abs() as f32
                / (average_surrounding_counts as f32);

            // NOTE: in the future, removing valid connections between nodes might not be desirable, as we might be able to resolve these
            // unitigs with e.g. larger k-value graphs.
            if rel_diff_0 > 0.2 && rel_diff_1 <= 0.2 {
                // 1 has similar counts as the neighbours
                // We remove the connections with 0
                // println!("here6");
                let tmpvec: Vec<EdgeIndex> = ptgraph
                    .graph
                    .edges_directed(midconns[0].0, Outgoing)
                    .chain(ptgraph.graph.edges_directed(midconns[0].0, Incoming))
                    .map(|er| er.id())
                    .collect();
                for e in tmpvec {
                    ptgraph.graph.remove_edge(e);
                }
            } else if rel_diff_0 <= 0.2 && rel_diff_1 > 0.2 {
                // 0 has similar counts as the neighbours
                // We remove the connections with 1
                // println!("here7");
                let tmpvec: Vec<EdgeIndex> = ptgraph
                    .graph
                    .edges_directed(midconns[1].0, Outgoing)
                    .chain(ptgraph.graph.edges_directed(midconns[1].0, Incoming))
                    .map(|er| er.id())
                    .collect();
                for e in tmpvec {
                    ptgraph.graph.remove_edge(e);
                }
            } else {
                // println!("here8");
                // Both have different counts as the neighbours
                // Perhaps both connections are unlikely. We remove all the connections of the bubble (i.e. we create a contig break).
                let tmpvec: Vec<EdgeIndex> = ptgraph
                    .graph
                    .edges_directed(midconns[0].0, Outgoing)
                    .chain(ptgraph.graph.edges_directed(midconns[0].0, Incoming))
                    .chain(
                        ptgraph
                            .graph
                            .edges_directed(midconns[1].0, Outgoing)
                            .chain(ptgraph.graph.edges_directed(midconns[1].0, Incoming)),
                    )
                    .map(|er| er.id())
                    .collect();
                for e in tmpvec {
                    ptgraph.graph.remove_edge(e);
                }
            }
            // In these cases we exit as there is nothing else to do
            return true;
        } else {
            // println!("here9");
            // The difference in the number of k-mers is very small, this could be a SNP. Let's pick up the highest count one.
            if node0w.counts > node1w.counts {
                // println!("here10");
                chosennode = 0;
                savedmidw = node0w.clone();
            } else if node0w.counts < node1w.counts {
                // println!("here11");
                chosennode = 1;
                savedmidw = node1w.clone();
            } else if node0w.abs_ind.len() > node1w.abs_ind.len() {
                // println!("here12");
                chosennode = 0;
                savedmidw = node0w.clone();
            } else {
                // println!("here13");
                chosennode = 1;
                savedmidw = node1w.clone();
            }
        }
    }

    // // Decision time! (simple)
    // if node0w.counts > node1w.counts {
    //     chosennode = 0;
    //     savedmidw = node0w.clone();
    // } else if node0w.counts < node1w.counts {
    //     chosennode = 1;
    //     savedmidw = node1w.clone();
    // } else if node0w.abs_ind.len() > node1w.abs_ind.len() {
    //     chosennode = 0;
    //     savedmidw = node0w.clone();
    // } else {
    //     chosennode = 1;
    //     savedmidw = node1w.clone();
    // }

    // We need to know if we are arriving, at the exit node of the bubble, to the Max ct.
    let midnodect = midconns[chosennode].1.get_from_and_to().1;
    let midconn2 = ptgraph.out_neighbours_bi(midconns[chosennode].0, midnodect)[0];
    let outct = midconn2.1.get_from_and_to().1;
    let outconn = ptgraph.out_neighbours_bi(midconn2.0, outct)[0];
    let savedoutw = ptgraph.graph.node_weight(midconn2.0).unwrap().clone();

    // By construction, the start node and the end node MUST have only one k-mer. If not, something is wrong.
    if ptgraph.graph.node_weight(startn).unwrap().abs_ind.len() > 1
        || ptgraph.graph.node_weight(midconn2.0).unwrap().abs_ind.len() > 1
    {
        panic!("Trying to remove a bubble with a start or end node with more than one k-mer!");
    };

    // We have copied the weights and have the information we need to continue, we can thus erase all the
    // nodes except from the start node.
    ptgraph.graph.remove_node(midconns[0].0); // intermediate node 0
    ptgraph.graph.remove_node(midconns[1].0); // intermediate node 1
    ptgraph.graph.remove_node(midconn2.0); // end node
                                           // At the end, the final node (the original start node) will have an internal node with and edge
                                           // MinToMin, as we are starting from the Min ct in any case. However, we might need to adjust the
                                           // following edges from the end node if we exit through the Max one.

    let mutrefw = ptgraph.graph.node_weight_mut(startn).unwrap();

    mutrefw.merge(&savedmidw, midconns[chosennode].1); // Intermediate node k-mer(s)
    mutrefw.set_internal_edge(EdgeType::MinToMin); // Inner edge
    mutrefw.abs_ind.push(savedoutw.abs_ind[0]); // Out node k-mer
    mutrefw.set_mean_counts(&[mutrefw.counts, savedmidw.counts, savedoutw.counts]);

    match outct {
        CarryType::Min => {
            // Now, we must link the start node with the following node to whatever node the end
            // node was connected.
            ptgraph
                .graph
                .add_edge(startn, outconn.0, EmptyEdge { t: outconn.1 });
            ptgraph
                .graph
                .add_edge(outconn.0, startn, EmptyEdge { t: outconn.1.rev() });
        }
        CarryType::Max => {
            // We must check what edges we had
            let tmptype = EdgeType::from_carrytypes(CarryType::Min, outconn.1.get_from_and_to().1);
            ptgraph
                .graph
                .add_edge(startn, outconn.0, EmptyEdge { t: tmptype });
            ptgraph
                .graph
                .add_edge(outconn.0, startn, EmptyEdge { t: tmptype.rev() });
        }
    }
    true
}

/// Remove dead input path.
#[inline]
fn remove_paths(ptgraph: &mut PtGraph, to_remove: Drain<NodeIndex>) {
    log::trace!("Removing {} dead paths", to_remove.len());
    for n in to_remove {
        ptgraph.graph.remove_node(n);
    }
}

/// Check if vertex initializes a dead path.
#[inline]
fn check_dead_path(
    ptgraph: &PtGraph,
    vertex: NodeIndex,
    output_vec: &mut Vec<NodeIndex>,
    k: usize,
    carryedge: EdgeType,
) {
    let mut current_vertex = vertex;
    output_vec.push(current_vertex);
    let cntopt = ptgraph.graph.node_weight(current_vertex);
    let mut cnt: usize;
    if let Some(thecntopt) = cntopt {
        cnt = thecntopt.abs_ind.len();
    } else {
        return;
    }

    let (mut ty, _) = carryedge.get_from_and_to();
    //     println!("> Checking tip from vertex {:?} with in_deg. {}, out_deg {}, and in the direction {:?}; k = {}",
    //         current_vertex, ptgraph.in_degree(current_vertex), ptgraph.out_degree(current_vertex),
    //         second_direction, k);
    let minnts = 100; // independent of this value, the minimum number of nts will be always k
    let limit = max(0, minnts - k + 1);
    // println!("\n START");
    loop {
        //         println!("\t- Iteration: curr_vertex {:?}, cnt {}", current_vertex, cnt);
        if cnt >= limit {
            // this path is not dead
            // println!("\t\t# Path deemed NOT dead! Count {}", cnt);
            output_vec.clear();
            return;
        }

        // This, when we change back to standard graphs surely will be more optimised. Now,
        // I don't know if we can do much better...
        let fwdneigh = ptgraph.out_neighbours_bi(current_vertex, ty); // NOTE: as all shrunk nodes have straight
                                                                      // internal edges, we don't need to check them,
                                                                      // because the entry and exit types will
                                                                      // always be the same.

        // if let Some(id) = ptgraph.graph.node_weight(current_vertex).unwrap().innerdir {
        //     println!("{:?}", id);
        // }
        let nfwdn = fwdneigh.len();

        if nfwdn != 1 {
            panic!("Not expected!");
        }

        let candidate_node = fwdneigh[0].0;
        let candidate_ty = fwdneigh[0].1.get_from_and_to().1;
        let bkgneigh_c = ptgraph.in_neighbours_bi(candidate_node, candidate_ty);
        let fwdneigh_c = ptgraph.out_neighbours_bi(candidate_node, candidate_ty);
        let nfwdn_c = fwdneigh_c.len();
        let nbkgn_c = bkgneigh_c.len();

        if nbkgn_c == 0 {
            panic!("Not expected! 2");
        } else if nbkgn_c != 1 {
            // We want to check before removing this tip that is the best one that could be removed.
            let mut altpath: Vec<Vec<NodeIndex>> = Vec::with_capacity(nbkgn_c - 1);
            let mut maxlen = 0;
            // println!("Last: {:?}, Vector: {:?}", *output_vec.last().unwrap(), bkgneigh_c);
            for n in bkgneigh_c.iter() {
                if n.0 == *output_vec.last().unwrap() {
                    // println!("This happens!");
                    continue;
                } else {
                    let mut tmppath: Vec<NodeIndex> = Vec::new();
                    check_backwards_path(
                        ptgraph,
                        n.0,
                        n.1.get_from_and_to().0,
                        &mut tmppath,
                        limit,
                    );
                    let tmplen = tmppath.len();
                    if tmplen != 0 {
                        altpath.push(tmppath);
                        if tmplen > maxlen {
                            maxlen = tmplen;
                        }
                    }
                }
            }

            // TODO/NOTE: perhaps add count critera instead of just length?
            // TODO: if the following if is not true, remove the correspondent tips and don't let that work for
            //       a different iteration.

            if maxlen != 0 && cnt > maxlen {
                // println!("Erasing other tip(s)! Count: {}, maxlen: {}", cnt, maxlen);
                // This means that we should NOT erase this tip, but the other(s), because they are shorter.
                output_vec.clear();
                for iv in altpath.iter_mut() {
                    output_vec.append(iv);
                }
            }
            return;
        }

        // Now we check the following neighbours
        if nfwdn_c == 0 {
            // Nowhere to continue after the candidate!
            // This is the situation (perfect sequence of kmers, no previous nor more forward neighbours in thec
            // candidate node) where we can just try to end this sequence of kmers and see whether we make the
            // minimum or not.

            output_vec.push(candidate_node);
            cnt += ptgraph
                .graph
                .node_weight(candidate_node)
                .unwrap()
                .abs_ind
                .len();

            if cnt >= limit {
                // this path is not dead
                // println!("\t\t# Path deemed NOT dead! Count {}", cnt);
                output_vec.clear();
            }
            return;
        } else if nfwdn_c == 1 {
            current_vertex = candidate_node;
            ty = candidate_ty;
            output_vec.push(current_vertex);
            cnt += ptgraph
                .graph
                .node_weight(current_vertex)
                .unwrap()
                .abs_ind
                .len();
        } else {
            // We don't want to remove anything in this case
            // println!("ALTERNATE ENDING");
            output_vec.clear();
            return;
        }
    }
}

fn check_backwards_path(
    ptgraph: &PtGraph,
    vertex: NodeIndex,
    mut ty: CarryType,
    output_vec: &mut Vec<NodeIndex>,
    kmerlimit: usize,
) {
    let mut current_vertex = vertex;
    output_vec.push(current_vertex);
    let mut cnt = ptgraph
        .graph
        .node_weight(current_vertex)
        .unwrap()
        .abs_ind
        .len();
    // println!("Starting search backwards from node {:?} and ty {:?}", vertex, ty);
    loop {
        if cnt >= kmerlimit {
            // this path is not dead
            output_vec.clear();
            return;
        }

        let bkgneigh = ptgraph.in_neighbours_bi(current_vertex, ty);
        let bkgneighlen = bkgneigh.len();
        if bkgneighlen == 0 {
            // Before just returning, let's check the outward neighbours:
            if ptgraph.out_neighbours_bi(current_vertex, ty).len() != 1 {
                // This is not another tip
                output_vec.clear();
            }
            return;
        } else if bkgneighlen == 1 {
            // BUG =================================
            // current_vertex = bkgneigh[0].0;
            // BUG =================================

            // Before confirming this new node, check if it is clearly member of a tip
            // if ptgraph.in_degree(current_vertex) > 1 || ptgraph.out_degree(current_vertex) > 1 {
            if ptgraph.out_neighbours_bi(current_vertex, ty).len() != 1 {
                // This is not another tip
                output_vec.clear();
                return;
            }

            current_vertex = bkgneigh[0].0;
            ty = bkgneigh[0].1.get_from_and_to().0;
            output_vec.push(current_vertex);
            cnt += ptgraph
                .graph
                .node_weight(current_vertex)
                .unwrap()
                .abs_ind
                .len();
        } else {
            // This is not another tip
            output_vec.clear();
            return;
        }
    }
}
