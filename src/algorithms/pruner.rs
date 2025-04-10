//! Various algorithms for graph pruning - removing unnecessary vertices/edges.

// use crate::graph_works::SEQUENCES;
use crate::graphs::pt_graph::{NodeIndex, PtGraph};
use crate::graphs::Graph;
use crate::EdgeWeight;
use crate::graphs::pt_graph::{CarryType, EdgeType};

use std::cmp::max;
// use std::process::exit;
use std::vec::Drain;


/// Describes prunable structure in the sense of genome assembly
pub trait Prunable: Clean {
    /// Remove all input and output dead paths
    fn remove_dead_paths(&mut self) -> bool;
    /// Find and remove all links that are impossible in bi-directed
    /// de Bruijn graphs derived from DNA sequences.
    fn remove_conflicting_links(&mut self) -> bool;
}

/// A trait for keeping the graph clean.
/// It keeps simple functions used for basic graph cleanups
pub trait Clean {
    /// Remove edges with weight below threshold.
    fn remove_weak_nodes(&mut self, threshold: EdgeWeight);
    /// Remove edges that are self-loops, i.e. those whose source
    /// and destination nodes are the same.
    fn remove_self_loops(&mut self);
}


impl Prunable for PtGraph {
    fn remove_dead_paths(&mut self) -> bool {
        log::info!("BEFORE PRUNING: {} nodes and {} edges",
            self.graph.node_count(),
            self.graph.edge_count());

        let mut dididoanything = false;
        log::info!("Starting graph pruning. Graph has {} externals, {} alone nodes, the remaining are internal.",
            self.externals_bi().len(),
            self.graph.node_indices().filter(|n| self.out_degree(*n) == 0 && self.in_degree(*n) == 0).count());

        let mut to_remove: Vec<NodeIndex> = vec![];
        loop {
            let mut path_check_vec = vec![];
            let externals : Vec<_> = self.externals_bi().into_iter().filter(|n| self.out_degree(*n) == 1).collect();

            // TODO: go back to "Graph" alone, for speed

            log::trace!("Detected {} externals", externals.len());

            for v in externals {
                check_dead_path(self,
                                   v,
                                   &mut path_check_vec,
                                   self.k,
                                   self.graph.edges_directed(v, petgraph::EdgeDirection::Outgoing).next().unwrap().weight().t);
                if !path_check_vec.is_empty() {
                    dididoanything = true;
                    to_remove.extend(path_check_vec.drain(..));
                }
            }

            // if there are no dead paths left pruning is done
            if to_remove.is_empty() {
                log::info!("Graph is pruned.");
                return dididoanything;
            }

            // reverse sort edge indices such that removal won't cause any troubles with swapped
            // edge indices (see `petgraph`'s explanation of `remove_edge`)
            // NOTE: this might be useful for the change to Graph backend, but now we don't need it.
            // to_remove.sort_by(|a, b| b.cmp(a));
            remove_paths(self, to_remove.drain(..));
        }
    }

    fn remove_conflicting_links(&mut self) -> bool {
        // Conflicting links: k-mers connected to others with a relatively larger count
        // TODO: apply

        return false;

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
    }
}


impl Clean for PtGraph {
    fn remove_weak_nodes(&mut self, threshold: EdgeWeight) {
        self.graph.retain_nodes(|g, n| g.node_weight(n).unwrap().counts >= threshold);
    }

    fn remove_self_loops(&mut self) {
        self.graph.retain_edges(|g, e| {let (n1, n2) = g.edge_endpoints(e).unwrap();
                                       n1 != n2});
    }
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
fn check_dead_path(ptgraph: &PtGraph, vertex: NodeIndex, output_vec: &mut Vec<NodeIndex>, k : usize, carryedge : EdgeType) {
    let mut current_vertex = vertex;
    output_vec.push(current_vertex);
    let cntopt = ptgraph.graph.node_weight(current_vertex);
    let mut cnt : usize;
    if cntopt.is_none() {
        return;
    } else {
        cnt = cntopt.unwrap().abs_ind.len();
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
        let fwdneigh;

        fwdneigh = ptgraph.out_neighbours_bi(current_vertex, ty); // NOTE: as all shrunk nodes have straight
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
        let candidate_ty   = fwdneigh[0].1.get_from_and_to().1;
        let bkgneigh_c = ptgraph.in_neighbours_bi(candidate_node, candidate_ty);
        let fwdneigh_c = ptgraph.out_neighbours_bi(candidate_node, candidate_ty);
        let nfwdn_c = fwdneigh_c.len();
        let nbkgn_c = bkgneigh_c.len();

        if nbkgn_c == 0 {
            panic!("Not expected! 2");
        } else if nbkgn_c != 1 {
            // We want to check before removing this tip that is the best one that could be removed.
            let mut altpath : Vec<Vec<NodeIndex>> = Vec::with_capacity(nbkgn_c - 1);
            let mut maxlen = 0;
            // println!("Last: {:?}, Vector: {:?}", *output_vec.last().unwrap(), bkgneigh_c);
            for n in bkgneigh_c.iter() {
                if n.0 == *output_vec.last().unwrap() {
                    // println!("This happens!");
                    continue;
                } else {
                    let mut tmppath : Vec<NodeIndex> = Vec::new();
                    check_backwards_path(&ptgraph, n.0, n.1.get_from_and_to().0, &mut tmppath, limit);
                    let tmplen = tmppath.len();
                    if tmplen != 0 {
                        altpath.push(tmppath);
                        if tmplen > maxlen {
                            maxlen  = tmplen;
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
            cnt += ptgraph.graph.node_weight(candidate_node).unwrap().abs_ind.len();

            if cnt >= limit {
                // this path is not dead
                // println!("\t\t# Path deemed NOT dead! Count {}", cnt);
                output_vec.clear();
            }
            return;

        } else if nfwdn_c == 1 {
            current_vertex  = candidate_node;
            ty              = candidate_ty;
            output_vec.push(current_vertex);
            cnt += ptgraph.graph.node_weight(current_vertex).unwrap().abs_ind.len();

        } else {
            // We don't want to remove anything in this case
            // println!("ALTERNATE ENDING");
            output_vec.clear();
            return;
        }
    }
}



fn check_backwards_path(ptgraph: &PtGraph, vertex: NodeIndex, mut ty: CarryType, output_vec: &mut Vec<NodeIndex>, kmerlimit : usize) {
    let mut current_vertex = vertex;
    output_vec.push(current_vertex);
    let mut cnt = ptgraph.graph.node_weight(current_vertex).unwrap().abs_ind.len();
    // println!("Starting search backwards from node {:?} and ty {:?}", vertex, ty);
    loop {
        if cnt >= kmerlimit {
            // this path is not dead
            output_vec.clear();
            return;
        }

        let bkgneigh    = ptgraph.in_neighbours_bi(current_vertex, ty);
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
            ty             = bkgneigh[0].1.get_from_and_to().0;
            output_vec.push(current_vertex);
            cnt += ptgraph.graph.node_weight(current_vertex).unwrap().abs_ind.len();
        } else {
            // This is not another tip
            output_vec.clear();
            return;
        }
    }
}
