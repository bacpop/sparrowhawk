//! Create string representation of contigs out of `Graph`.

use super::shrinker::Shrinkable;
use crate::graphs::Graph;
use super::pruner::Clean;
use crate::graphs::pt_graph::{CarryType, NodeIndex, NodeStruct, PtGraph};
use std::cmp::max;
// use std::process::exit;

use rayon::prelude::*;
use petgraph::EdgeDirection;
use petgraph::algo::{connected_components, tarjan_scc};
use petgraph::visit::EdgeRef;
use petgraph;


/// Collapse `Graph` into `SerializedContigs`.
pub trait Collapsable: Shrinkable {
    /// Collapses `Graph` into `SerializedContigs`.
    fn collapse(self, path : Option<&String>) -> SerializedContigs;
}


/// Representation of serialized contig.
pub type SerializedContig = Vec<NodeStruct>;


/// Collection of serialized contigs.
pub type SerializedContigs = Vec<Vec<NodeStruct>>;


impl Collapsable for PtGraph {
    fn collapse(mut self, path_ : Option<&String>) -> SerializedContigs {
        let mut contigs: SerializedContigs = vec![];

        log::info!("Removing self-loops (temporal restriction)");
        self.remove_self_loops();

        log::info!("Graph has {} weakly connected component(s), among which {} are single nodes.",
            connected_components(&petgraph::graph::Graph::from(self.graph.clone())),
            self.graph.node_indices().filter(|n| self.graph.neighbors_directed(*n, EdgeDirection::Outgoing).count() == 0 && self.graph.neighbors_directed(*n, EdgeDirection::Incoming).count() == 0).count()
        );

        log::info!("Starting collapse loop.");
        let minnts = 100; // independent of this value, the minimum number of nts will be always k
        let limit = max(0, minnts - self.k + 1);

        loop {
            loop {
                // get all starting nodes, i.e. nodes with in_degree == 0
                let externals = self.externals_bi();
                println!("\t- Loop over {} external nodes.", externals.len());
                if externals.len() == 0 {
                    break;
                }
                // create contigs from each starting node
                for n in externals {
                    // We need first to take care of perfect contigs, almost-already provided as such. These
                    // are seen as nodes with no incoming/outcoming edges.
                    if self.graph.contains_node(n) {
                        if self.get_good_connections_degree(n) == 0 {
                                println!("\t\t# Isolated node.");
                            let thecont = vec![self.graph.node_weight(n).unwrap().clone()];
                            if get_contig_length(&thecont) > limit {
                                contigs.push(thecont);
                            }
                            // contigs.push(thecont);
                            self.graph.remove_node(n);
        //                         tobreak = false;
                        } else {
                                // println!("\t\t# Non-isolated node: creating contig.");
                            stacker::maybe_grow(32 * 1024, 1 * 1024 * 1024, || {
                                let contigs_ = contigs_from_vertex(&mut self, n);
                                contigs.extend(contigs_.into_iter().filter(|c| get_contig_length(c) > limit).collect::<Vec<_>>());
                            });
                        }
                    }
                }
            }
//                 if tobreak {break};
            println!("\t- Loops over external nodes finished.");

            // break;

            let tmpc = self.graph.node_count();
            if tmpc != 0 { // we guarantee that there's at least one node to unwrap here
                println!("\t\t# {} nodes remain. Starting to build from middle node.", tmpc);
    //                 break;
    //                 println!("\t\t# Getting strongly-connected components");

                // log::info!("Saving remaining graph with loops/bubbles/circumferences/whatever as DOT file...");
                // if path_.is_some() && dosave {
                //     let mut wbuf = set_ostream(&Some(path_.unwrap().clone().replace(".dot", "_beforeattackingtheremains.dot")));
                //     self.write_to_dot(&mut wbuf);
                //     dosave = false;
                // }

                // This call to stacker::grow here is needed because of the algorithm that is run to obtain the
                // strongly-connected components. It is recursive, so in very entangled graphs (and/or when k is
                // low, i.e. k ~< 15), it might lead to a stack overflow.
                stacker::grow(100 * 1024 * 1024, || {
                    let sccvec : Vec<Vec<NodeIndex>> = tarjan_scc(&self.graph);
                    let node_in_cycle = sccvec[0].last().unwrap();

                    println!("\t\t# Remaining nodes {}, remaining SCCs {}, starting with {} neighbours",
                        tmpc,
                        sccvec.len(),
                        self.get_good_connections_degree(*node_in_cycle));

                    let thecontigs = contigs_from_intermediate_vertex(&mut self, *node_in_cycle);

                    contigs.extend(thecontigs.into_iter().filter(|c| get_contig_length(c) > limit).collect::<Vec<_>>());
                });
                println!("\t\t# Finished creating one contig from starting circle.");

            } else {
                break;
            }
        }

        log::trace!("{} nodes left in the graph after collapse", self.graph.node_count());
        log::info!("Collapse ended. Created {} contigs", contigs.len());

        contigs
    }
}


fn get_contig_length(vec : &Vec<NodeStruct>) -> usize {
    let mut outnum = 0;
    for ns in vec.iter() {
        outnum += ns.abs_ind.len();
    }
    outnum
}



// Main collapse function/method
#[inline]
fn contigs_from_vertex(ptgraph: &mut PtGraph, v: NodeIndex) -> SerializedContigs {
    let mut contigs: SerializedContigs = vec![];
    let mut contig : SerializedContig  = vec![];
    let mut current_vertex = v;
    // let mut current_edge_index;
    let mut target;
    let outeds : Vec<_>  = ptgraph.graph.edges_directed(v, EdgeDirection::Outgoing).map(|e| e.id()).collect();
    let mut current_type = ptgraph.graph.edge_weight(outeds[0]).unwrap().t.get_from_and_to().0;
    let mut outneighs    = ptgraph.out_neighbours_bi(v, current_type);
    let mut num_following = outneighs.len();
    let mut num_preceding = 0;


    // println!("\t\t\t% START! From vertex {:?}", current_vertex);
    loop {
        // println!("\t\t\t% Iteration start: vertex {:?} num_preceding {} num_following {}", current_vertex, num_preceding, num_following);

        if num_following == 1 && num_preceding == 0 {
        // Ok, so we can continue, let's go!
            // current_edge_index = outneighs[0].1;
            // println!("\t # Cont!");
        } else if num_following == 0 && num_preceding == 0 { // We're finishing!!
//             println!("{:?}", current_vertex);

            let mut nwtocopy = ptgraph.graph.node_weight(current_vertex).unwrap().clone();
            if let Some(innvtx) = nwtocopy.innerdir {
                if current_type != innvtx.get_from_and_to().0 {
                    nwtocopy.abs_ind.reverse();
                }
            }

            contig.push(nwtocopy);
            // contig.push(ptgraph.graph.node_weight(current_vertex).expect(&format!("Not found vertex {:?}, that should have num_following {} and num_preceding {}", current_vertex, num_following, num_preceding)).clone());

//             println!("We are stopping. This is the last contig: {:?}", contig);
            contigs.push(contig.clone());
            contig.clear();
//             decrease_weight(ptgraph, current_vertex);
            ptgraph.graph.remove_node(current_vertex);
            // println!("\t # Ending without possible continuation!");
            return contigs;
        } else {
            // We've found an ambiguous node/bifurcation, thus we need to stop the current contig and clear the vector
            contigs.push(contig.clone());
            contig.clear();

            // And now what we do depends on the neighbours from this new vertex. OR NOT: LET'S FINISH FOR NOW!
            if num_following == 0 {
                // We cannot continue.
                // println!("\t # Ending because we reach an ambiguous node that is a sink!");

                // TEST BEGIN
                ptgraph.graph.remove_node(current_vertex);
                // TEST END

                return contigs;
            }


            // TEST BEGIN
            ptgraph.graph.remove_node(current_vertex);
            return contigs;
            // TEST END
        }


        // If we arrived here, current_vertex is either considered good to be added to the current
        // contig, or we have created a contig break and we are starting from this ambiguous node
        // and also current_edge_index is the vertex through which we should continue our
        // journey, or we have either a simple loop or a circumference to deal with

        // We add the current_vertex to the contig
        let mut nwtocopy = ptgraph.graph.node_weight(current_vertex).unwrap().clone();
        if let Some(innvtx) = nwtocopy.innerdir {
            if current_type != innvtx.get_from_and_to().0 {
                nwtocopy.abs_ind.reverse();
            }
        }

        contig.push(nwtocopy);

        // We get our next soon-to-be current_vertex (now named "target")
        // (_, target) = ptgraph.graph.edge_endpoints(current_edge_index).unwrap();
        target = outneighs[0].0;

//         println!("\t\t\t% Current vertex: {:?}, Target: {:?}", current_vertex, target);
        if current_vertex == target {panic!("FATAL: continuing to the same vertex!")};


        // We update the variables and get ready to do another iteration!
        current_type = outneighs[0].1.get_from_and_to().1;
        ptgraph.graph.remove_node(current_vertex);
        current_vertex = target;
        num_preceding = ptgraph.in_degree_bi(current_vertex, current_type);
        outneighs = ptgraph.out_neighbours_bi(current_vertex, current_type);
        num_following = outneighs.len();
//         println!("{} {} {:?}", num_preceding, num_following, target);
    }
}


#[inline]
fn contigs_from_intermediate_vertex(ptgraph: &mut PtGraph, v: NodeIndex) -> SerializedContigs {
    let mut contigs: SerializedContigs = vec![];
    let mut contig : SerializedContig  = vec![];
    let mut current_vertex = v;
    // let mut current_edge_index;
    let mut target;

    // We need to get the carrytype, the edges, and so on before we can begin. We'll try to set them to get a forward
    // direction with only one neighbour, if possible.
    let outeds;
    let mut outmin = Vec::new();
    let mut outmax = Vec::new();
    // let mut started = false;

    for e in ptgraph.graph.edges_directed(v, EdgeDirection::Outgoing) {
        if e.weight().t.get_from_and_to().0 == CarryType::Min {
            outmin.push(e.id());
        } else {
            outmax.push(e.id());
        }
    }
    let outminlen = outmin.len();
    let outmaxlen = outmax.len();
    match (outminlen, outmaxlen) {
        (0, 0) | (0, 1) | (1, 0) => panic!("External node!!!"),
        (1, 1) | (1, _) => outeds = outmin, // We select the minimum outgoing edges
        (_, 1)          => outeds = outmax, // We select the maximum outgoing edges
        (_, _)          => { // We check whether they are the same and, if not, we select the first id from the minimum (this is clearly improvable)
            if outminlen <= outmaxlen {
                outeds = outmin;
            } else {
                outeds = outmax;
            }
        },
    }

    let mut current_type    = ptgraph.graph.edge_weight(outeds[0]).unwrap().t.get_from_and_to().0;
    let mut outneighs       = ptgraph.out_neighbours_bi(v, current_type);
    let mut num_following   = outneighs.len();
    let mut num_preceding = 0; /////// This is strictly speaking always false here, but it is only for the first iteration.
                                    // Afterwards, we respect its true value to decide whether we stop or not the contig formation.

    // println!("\t\t\t% START! From vertex {:?}", current_vertex);
    loop {
        // println!("\t\t\t% Iteration start: vertex {:?} num_preceding {} num_following {}", current_vertex, num_preceding, num_following);

        if num_following == 1 && num_preceding == 0 {
        // Ok, so we can continue, let's go!
            // current_edge_index = outneighs[0].1;
            // println!("\t # Cont!");
        } else if num_following == 0 && num_preceding == 0 { // We're finishing!!
//             println!("{:?}", current_vertex);

            let mut nwtocopy = ptgraph.graph.node_weight(current_vertex).unwrap().clone();
            if let Some(innvtx) = nwtocopy.innerdir {
                if current_type != innvtx.get_from_and_to().0 {
                    nwtocopy.abs_ind.reverse();
                }
            }

            contig.push(nwtocopy);
            // contig.push(ptgraph.graph.node_weight(current_vertex).expect(&format!("Not found vertex {:?}, that should have num_following {} and num_preceding {}", current_vertex, num_following, num_preceding)).clone());

//             println!("We are stopping. This is the last contig: {:?}", contig);
            contigs.push(contig.clone());
            contig.clear();
//             decrease_weight(ptgraph, current_vertex);
            ptgraph.graph.remove_node(current_vertex);
            // println!("\t # Ending without possible continuation!");
            return contigs;
        } else {
            // We've found an ambiguous node/bifurcation, thus we need to stop the current contig and clear the vector
            contigs.push(contig.clone());
            contig.clear();

            // And now what we do depends on the neighbours from this new vertex. OR NOT: LET'S FINISH FOR NOW!
            if num_following == 0 {
                // We cannot continue.
                // println!("\t # Ending because we reach an ambiguous node that is a sink!");
                return contigs;
            }
        }


        // If we arrived here, current_vertex is either considered good to be added to the current
        // contig, or we have created a contig break and we are starting from this ambiguous node
        // and also current_edge_index is the vertex through which we should continue our
        // journey, or we have either a simple loop or a circumference to deal with

        // We add the current_vertex to the contig
        let mut nwtocopy = ptgraph.graph.node_weight(current_vertex).unwrap().clone();
        if let Some(innvtx) = nwtocopy.innerdir {
            if current_type != innvtx.get_from_and_to().0 {
                nwtocopy.abs_ind.reverse();
            }
        }

        contig.push(nwtocopy);

        // We get our next soon-to-be current_vertex (now named "target")
        // (_, target) = ptgraph.graph.edge_endpoints(current_edge_index).unwrap();
        target = outneighs[0].0;

//         println!("\t\t\t% Current vertex: {:?}, Target: {:?}", current_vertex, target);
        if current_vertex == target {panic!("FATAL: continuing to the same vertex!")};


        // We update the variables and get ready to do another iteration!
        current_type = outneighs[0].1.get_from_and_to().1;
        ptgraph.graph.remove_node(current_vertex);
        current_vertex = target;
        num_preceding = ptgraph.in_degree_bi(current_vertex, current_type);
        outneighs = ptgraph.out_neighbours_bi(current_vertex, current_type);
        num_following = outneighs.len();
//         println!("{} {} {:?}", num_preceding, num_following, target);
    }

}
