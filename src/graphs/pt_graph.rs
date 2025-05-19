use crate::algorithms::builder::{Build, Init};
use super::Graph;

use crate::Idx;

use core::fmt;
use std::io::Write;
use std::collections::BTreeSet;


use petgraph;
use petgraph::visit::EdgeRef;
use petgraph::dot::{Dot, Config};

use nohash_hasher::NoHashHasher;
use petgraph::Direction::Outgoing;
use std::{collections::HashMap, hash::BuildHasherDefault, cell::*};
use crate::HashInfoSimple;

/// Type denoting index of edge.
pub type EdgeIndex = petgraph::stable_graph::EdgeIndex<Idx>;
/// Type denoting index of node.
pub type NodeIndex = petgraph::stable_graph::NodeIndex<Idx>;


/// Class of a petgraph-based graph in which we'll implement the graphs traits defined in other
/// files of the crate.
#[derive(Default)]
pub struct PtGraph {
    /// Petgraph stable-graph: foundations of the graph structure.
    pub graph : petgraph::stable_graph::StableGraph<NodeStruct, EmptyEdge, petgraph::Directed, Idx>,

    /// k-value of the graph.
    pub k     : usize,
}


/// Structure that contains the information of a node.
#[derive(Clone,Debug)]
pub struct NodeStruct {
    /// Value that reflects either the counts of one k-mer, or a
    /// proxy value for shrunk nodes.
    pub counts   : u16,

    /// List of hashes os k-mers
    pub abs_ind  : Vec<u64>,

    /// Inner edge for those edges result of a shrinkage
    pub innerdir : Option<EdgeType>,
}


impl NodeStruct {
    /// Allows merging two nodes and their inner information
    pub fn merge(&mut self, other: &NodeStruct, tytoother : EdgeType) {   // NEW
        // Note that counts still need to be recalculated BEFORE!!!!

        // First, we'll see if we have an internal edge already.
        if let Some(thisid) = self.innerdir {
            // This is a bit more complicated...
            let thisct = thisid.get_from_and_to().0;

            if thisct == tytoother.get_from_and_to().0 {
                // Great! We don't need to change nothing from this vector.
                // Let's see now the OTHER vector...
                if let Some(otherid) = other.innerdir {
                    // This is a bit more complicated...
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        // We can directly merge the vectors
                        self.abs_ind.extend(&other.abs_ind);
                    } else {
                        // This is a bit more complicated: we need to INVERT
                        // the other vector BEFORE joining them
                        let mut newvec = other.abs_ind.clone();
                        newvec.reverse();
                        self.abs_ind.extend(&newvec);
                    }
                } else {
                    // This is easy!
                    self.abs_ind.extend(&other.abs_ind);
                }

            } else {
                // This is not great, this node is currently in the opposite orientation to that
                // of the current merge (as defined by tytoother).
                // Thus, we need to insert whatever comes at the beginning of the vector

                if let Some(otherid) = other.innerdir {
                    // This is a bit more complicated... // TODO
                    let otherct = otherid.get_from_and_to().0;
                    if tytoother.get_from_and_to().1 == otherct {
                        // The other shrunk node is aligned with the edge that connects it with this one.
                        // As ours is not, we must prepend the REVERSED information of the other shrunk node.
                        let mut tmpvec = other.abs_ind.clone();
                        tmpvec.reverse();
                        self.abs_ind.splice(0..0, tmpvec.iter().cloned());
                    } else {
                        // This is a bit more complicated: we DON'T need to reverse
                        // the other vector BEFORE joining them, because they are already in the same order.
                        self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                    }
                } else {
                    self.abs_ind.splice(0..0, other.abs_ind.iter().cloned());
                }
            }

        } else {
            // This is easier!
            // Let's see if the other node has an internal edge itself

            if let Some(otherid) = other.innerdir {
                // This is a bit more complicated...
                let otherct = otherid.get_from_and_to().0;
                if tytoother.get_from_and_to().1 == otherct {
                    // We can directly merge the vectors
                    self.abs_ind.extend(&other.abs_ind);
                } else {
                    // This is a bit more complicated: we need to INVERT
                    // the other vector BEFORE joining them
                    let mut newvec = other.abs_ind.clone();
                    newvec.reverse();
                    self.abs_ind.extend(&newvec);
                }
            } else {
                // This is very easy! We just have to join the vectors
                self.abs_ind.extend(&other.abs_ind);
            }
        }
    }


    /// Checks whether it is needed to reverse the list of hashes and the inner edge.
    pub fn invert_if_needed(&mut self, outedge : EdgeType) {
        if let Some(id) = self.innerdir {
            if id.get_from_and_to().1 != outedge.get_from_and_to().0 {
                self.abs_ind.reverse();
                self.innerdir.unwrap().rev();
            }
        }
    }


    /// Sets the counts of the node as the mean of the vector you give the function.
    pub fn set_mean_counts(&mut self, countsvec: &Vec<u16>) {
        self.counts = (countsvec.iter().map(|&e| e as u32).sum::<u32>() as f32 / countsvec.len() as f32).round() as u16;
        // Note that counts still need to be recalculated BEFORE!!!!
    }


    /// Sets the type of the internal edge as the one you provide the function.
    pub fn set_internal_edge(&mut self, ed : EdgeType) {
        // self.innerdir = Some(ed);
        match ed {
            EdgeType::MinToMin | EdgeType::MaxToMax => self.innerdir = Some(ed),
            _ => panic!("Non-valid internal edge!"),
        }
    }
}


impl fmt::Display for NodeStruct {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.counts)
    }
}



/// Enum that describes the type of one edge of the graph (essentially,
/// from which hash it comes (either canonical/minimum or non-canonical/maximum)
/// and with what it is linked (again, either min/max).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeType {
    /// Links canonical hash to canonical hash
    MinToMin,
    /// Links non-canonical hash to non-canonical hash
    MaxToMax,
    /// Links canonical hash to non-canonical hash
    MinToMax,
    /// Links non-canonical hash to canonical hash
    MaxToMin,
}


impl EdgeType {
    /// Reverses the EdgeType, changing it to the DNA de Bruijn-graph edge that would exist
    /// and begin at the end of the original EdgeType and end at the beginning of the same one.
    pub fn rev(&self) -> EdgeType {
        match self {
            EdgeType::MinToMin => EdgeType::MaxToMax,
            EdgeType::MaxToMax => EdgeType::MinToMin,
            EdgeType::MinToMax | EdgeType::MaxToMin => *self,
        }
    }

    /// Returns the CarryType of the origin and end of the edges, i.e. the canonicality of the
    /// hashes that this edge connects.
    pub fn get_from_and_to(&self) -> (CarryType, CarryType) {
        match self {
            EdgeType::MinToMin => (CarryType::Min, CarryType::Min),
            EdgeType::MaxToMax => (CarryType::Max, CarryType::Max),
            EdgeType::MaxToMin => (CarryType::Max, CarryType::Min),
            EdgeType::MinToMax => (CarryType::Min, CarryType::Max),
        }
    }

    /// Constructs an EdgeType from the canonicality (CarryType) of the source and end of the edges.
    pub fn from_carrytypes(first : CarryType, second :  CarryType) -> EdgeType {
        match (first, second) {
            (CarryType::Min, CarryType::Min) => EdgeType::MinToMin,
            (CarryType::Max, CarryType::Max) => EdgeType::MaxToMax,
            (CarryType::Min, CarryType::Max) => EdgeType::MinToMax,
            (CarryType::Max, CarryType::Min) => EdgeType::MaxToMin,
        }
    }

    /// Answers whether the EdgeType is of the direct types (i.e. linkes a canonical hash with another canonical one,
    /// or the equivalent with non-canonicals).
    pub fn is_direct(&self) -> bool {
        match self {
            EdgeType::MinToMin | EdgeType::MaxToMax => true,
            _                                       => false,
        }
    }
}


/// This enum describes the canonicality of a hash
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarryType {
    /// Canonical hash (the minimum hash of the pair)
    Min,
    /// Non-canonical hash (the maximum hash of the pair)
    Max,
}


/// This struct contains the information of an edge in the graph, which is "empty" because it only contains the type of the edge
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EmptyEdge {
    /// Type of the edge.
    pub t : EdgeType,
}


impl fmt::Display for EmptyEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "")
    }
}


impl Graph for PtGraph {
    type NodeIdentifier = NodeIndex;
    type AmbiguousNodes = BTreeSet<NodeIndex>;

    #[inline]
    fn get_ambiguous_nodes(&self) -> Self::AmbiguousNodes {
        self.graph.node_indices()
            .filter(|n| {
                let in_degree = self.in_degree(*n);
                let out_degree = self.out_degree(*n);
                in_degree > 1 || out_degree > 1 || (in_degree == 0 && out_degree == 1)
            })
            .collect::<Self::AmbiguousNodes>()
    }

    #[inline]
    fn get_ambiguous_nodes_bi(&self) -> Self::AmbiguousNodes {
        self.graph.node_indices()
            .filter(|n| {
                let conns = self.get_good_connections_degree(*n);

                // TODO: this conditional could be merged easily in the get_good_connections_degree call,
                // making everything more efficient.

                if conns == 0 {
                    return false;
                } else if conns == 2 {
                    // If conns == 2, we just need to avoid a perfect intermediate kmer in a sequence of them.
                    if self.out_degree_min(*n) == 1 {
                        return false;
                    }
                }

                return true;
            })
            .collect::<Self::AmbiguousNodes>()
    }


    #[inline]
    fn get_good_neighbours_bi(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)> {
        self.graph.edges_directed(node, Outgoing).filter(|e| e.source() != e.target()).map(|e| (e.target(), e.weight().t)).collect()
    }


    #[inline]
    fn get_good_connections_degree(&self, node: Self::NodeIdentifier) -> usize {
        // TODO: optimise this, this most surely is very slow.
        // NOTE: this does not consider self-loops.
        self.graph.edges_directed(node, Outgoing).filter(|e| e.source() != e.target()).map(|e| (e.target(), e.weight().t)).count()
    }

    #[inline]
    fn node_has_self_loops(&self, node: Self::NodeIdentifier) -> bool {
        for e in self.graph.edges_directed(node, Outgoing) {
            if e.target() == node {
                return true;
            }
        }

        // NOTE: add also the incoming???

        return false;
    }

    #[inline]
    fn out_degree(&self, node: Self::NodeIdentifier) -> usize {
        self.graph.neighbors_directed(node, petgraph::EdgeDirection::Outgoing)
            .count()
    }

    #[inline]
    fn in_degree_bi(&self, node: Self::NodeIdentifier, ty :CarryType) -> usize {
        match ty {
            CarryType::Min => self.in_degree_min(node),
            CarryType::Max => self.in_degree_max(node),
        }
    }

    #[inline]
    fn out_degree_bi(&self, node: Self::NodeIdentifier, ty :CarryType) -> usize {
        match ty {
            CarryType::Min => self.out_degree_min(node),
            CarryType::Max => self.out_degree_max(node),
        }
    }

    #[inline]
    fn out_degree_min(&self, node: Self::NodeIdentifier) -> usize {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Outgoing).collect::<Vec<_>>();
        let mut count = 0;
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MinToMin | EdgeType::MinToMax => count += 1,
                _ => (),
            }
        }

        count
    }

    #[inline]
    fn out_degree_max(&self, node: Self::NodeIdentifier) -> usize {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Outgoing).collect::<Vec<_>>();
        let mut count = 0;
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MaxToMax | EdgeType::MaxToMin => count += 1,
                _ => (),
            }
        }

        count
    }

    #[inline]
    fn in_neighbours_bi(&self, node: Self::NodeIdentifier, ty : CarryType) -> Vec<(NodeIndex, EdgeType)> {
        match ty {
            CarryType::Min => self.in_neighbours_min(node),
            CarryType::Max => self.in_neighbours_max(node),
        }
    }

    #[inline]
    fn out_neighbours_bi(&self, node: Self::NodeIdentifier, ty : CarryType) -> Vec<(NodeIndex, EdgeType)> {
        match ty {
            CarryType::Min => self.out_neighbours_min(node),
            CarryType::Max => self.out_neighbours_max(node),
        }
    }

    #[inline]
    fn out_neighbours_min(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)> {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Outgoing).collect::<Vec<_>>();
        let mut outneigh =Vec::new();
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MinToMin | EdgeType::MinToMax => outneigh.push( (ie.target(), ie.weight().t) ),
                _ => (),
            }
        }
        outneigh
    }

    #[inline]
    fn out_neighbours_max(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)> {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Outgoing).collect::<Vec<_>>();
        let mut outneigh =Vec::new();
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MaxToMin | EdgeType::MaxToMax => outneigh.push( (ie.target(), ie.weight().t) ),
                _ => (),
            }
        }
        outneigh
    }

    #[inline]
    fn in_degree(&self, node: Self::NodeIdentifier) -> usize {
        self.graph.neighbors_directed(node, petgraph::EdgeDirection::Incoming)
            .count()
    }

    #[inline]
    fn in_degree_min(&self, node: Self::NodeIdentifier) -> usize {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Incoming).collect::<Vec<_>>();
        let mut count = 0;
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MaxToMin | EdgeType::MinToMin => count += 1,
                _ => (),
            }
        }
        count
    }

    #[inline]
    fn in_degree_max(&self, node: Self::NodeIdentifier) -> usize {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Incoming).collect::<Vec<_>>();
        let mut count = 0;
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MinToMax | EdgeType::MaxToMax => count += 1,
                _ => (),
            }
        }
        count
    }

    #[inline]
    fn in_neighbours_min(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)> {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Incoming).collect::<Vec<_>>();
        let mut outneigh =Vec::new();
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MinToMin | EdgeType::MaxToMin => outneigh.push( (ie.source(), ie.weight().t) ),
                _ => (),
            }
        }
        outneigh
    }

    #[inline]
    fn in_neighbours_max(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)> {
        let tmpoutneigh = self.graph.edges_directed(node, petgraph::EdgeDirection::Incoming).collect::<Vec<_>>();
        let mut outneigh =Vec::new();
        for ie in tmpoutneigh.iter() {
            match ie.weight().t {
                EdgeType::MinToMax | EdgeType::MaxToMax => outneigh.push( (ie.source(), ie.weight().t) ),
                _ => (),
            }
        }
        outneigh
    }

    #[inline]
    fn externals_bi(&self) -> Vec<NodeIndex> {
        self.graph.node_indices().filter(|n| {
            // self.out_degree(*n) == 1 && self.in_degree(*n) == 1
            let mut it = self.graph.edges_directed(*n, petgraph::Direction::Incoming);
            let ct : CarryType;
            if let Some(e) = it.next() {
                ct = e.weight().t.get_from_and_to().1;
            } else {

                // println!("jojojojo");
                return true;
            }

            while let Some(e) = it.next() {
                if e.weight().t.get_from_and_to().1 != ct {
                    return false;
                }
            }
            return true;
        }).collect::<Vec<NodeIndex>>()
    }

    #[inline]
    fn externals_bi_min(&self, thedir : petgraph::EdgeDirection) -> Vec<NodeIndex> { // Makes no sense, REMOVE

        self.graph.node_indices().filter(|n| {
            if self.out_degree(*n) != 1 || self.in_degree(*n) != 1 {
                return false
            }
            match self.graph.edges_directed(*n, petgraph::EdgeDirection::Outgoing).next().unwrap().weight().t {
                EdgeType::MinToMin => match thedir {
                    petgraph::EdgeDirection::Outgoing => false,
                    petgraph::EdgeDirection::Incoming => true,
                },
                EdgeType::MinToMax => match thedir {
                    petgraph::EdgeDirection::Outgoing => false,
                    petgraph::EdgeDirection::Incoming => true,
                },
                EdgeType::MaxToMin => match thedir {
                    petgraph::EdgeDirection::Outgoing => true,
                    petgraph::EdgeDirection::Incoming => false,
                },
                EdgeType::MaxToMax => match thedir {
                    petgraph::EdgeDirection::Outgoing => true,
                    petgraph::EdgeDirection::Incoming => false,
                },
            }
        }).collect::<Vec<NodeIndex>>()
    }

    #[inline]
    fn externals_bi_max(&self, thedir : petgraph::EdgeDirection) -> Vec<NodeIndex> { // Makes no sense, REMOVE
        self.graph.node_indices().filter(|n| {
            if self.out_degree(*n) != 1 || self.in_degree(*n) != 1 {
                return false
            }
            match self.graph.edges_directed(*n, petgraph::EdgeDirection::Outgoing).next().unwrap().weight().t {
                EdgeType::MinToMin => match thedir {
                    petgraph::EdgeDirection::Outgoing => true,
                    petgraph::EdgeDirection::Incoming => false,
                },
                EdgeType::MinToMax => match thedir {
                    petgraph::EdgeDirection::Outgoing => true,
                    petgraph::EdgeDirection::Incoming => false,
                },
                EdgeType::MaxToMin => match thedir {
                    petgraph::EdgeDirection::Outgoing => false,
                    petgraph::EdgeDirection::Incoming => true,
                },
                EdgeType::MaxToMax => match thedir {
                    petgraph::EdgeDirection::Outgoing => false,
                    petgraph::EdgeDirection::Incoming => true,
                },
            }
        }).collect::<Vec<NodeIndex>>()
    }

    fn write_to_dot<W: Write>(&self, f: &mut W) {
//         println!("{}", Dot::new(&self.graph));
//         fs::write(path_, Dot::new(&self.graph)).expect("Unable to write file");
        let mut graphfordot = self.graph.clone();
        graphfordot.retain_nodes(|g, n| g.neighbors_directed(n, petgraph::EdgeDirection::Outgoing).count() != 0 || g.neighbors_directed(n, petgraph::EdgeDirection::Incoming).count() != 0);
        let output = format!("{:?}", Dot::with_attr_getters(&graphfordot,
            &[Config::NodeNoLabel, Config::EdgeNoLabel],
            // &|_, e| format!("label = \"{:?}\"", e.weight().t),
            // &|_, n| format!("label = \"{:?} | {:?}\"", n.1.counts, n.1.abs_ind.len()),));
            &|_, e| format!("label = \"{:?}\"", e.weight().t),
            &|_, n| format!("label = \"{:?} | {:?}\" counts = {:?} kmers = {:?}", n.1.counts, n.1.abs_ind.len(), n.1.counts, n.1.abs_ind.len()),));

        let _ = f.write(&output.into_bytes()[..]);
    }

    fn get_dot_string(&self) -> String {
//         println!("{}", Dot::new(&self.graph));
//         fs::write(path_, Dot::new(&self.graph)).expect("Unable to write file");
        let mut graphfordot = self.graph.clone();
        graphfordot.retain_nodes(|g, n| g.neighbors_directed(n, petgraph::EdgeDirection::Outgoing).count() != 0 || g.neighbors_directed(n, petgraph::EdgeDirection::Incoming).count() != 0);
        let output = format!("{:?}", Dot::with_attr_getters(&graphfordot,
            &[Config::NodeNoLabel, Config::EdgeNoLabel],
            // &|_, e| format!("label = \"{:?}\"", e.weight().t),
            // &|_, n| format!("label = \"{:?} | {:?}\"", n.1.counts, n.1.abs_ind.len()),));
            &|_, e| format!("label = \"{:?}\"", e.weight().t),
            &|_, n| format!("label = \"{:?} | {:?}\" counts = {:?} kmers = {:?}", n.1.counts, n.1.abs_ind.len(), n.1.counts, n.1.abs_ind.len()),));

        return output;
    }

    fn write_to_gfa<W: Write>(&self, f: &mut W) {
        let mut output : String = "VN:Z:2.0\n".to_owned();

        // Given the duality of the graph, to do the exportation of it as GFA, it is enough to assume that e.g. the canonical hashes represent the direct strand.
        // With that, the resulting nodes and edges will represent, by construction, correctly both strands.

        // Add nodes
        self.graph.node_indices().map(|ni| {
            output.push_str( format!("S\t{}\t{}\t*\n", ni.index(), self.k + self.graph.node_weight(ni).unwrap().abs_ind.len() - 1).as_str() );
        } );

        // Add edges
        self.graph.edge_indices().map(|ei| {
            let (sid, tid) = self.graph.edge_endpoints(ei).unwrap();
            let (st, tt)   = self.graph.edge_weight(ei).unwrap().t.get_from_and_to();

            let ssign : &str;
            let tsign : &str;
            let sbeg  : String;
            let send  : String;
            let tbeg  : String;
            let tend  : String;

            match st {
                CarryType::Min => {
                    ssign = "+";

                    let tmplen = self.graph.node_weight(sid).unwrap().abs_ind.len();
                    sbeg = format!("{}",  tmplen);
                    send = format!("{}$", self.k + tmplen - 1);
                },
                CarryType::Max => {
                    ssign = "-";

                    sbeg = format!("{}", 0);
                    send = format!("{}", self.k - 1);
                },
            }

            match tt {
                CarryType::Min => {
                    tsign = "+";
                    tbeg = format!("{}", 0);
                    tend = format!("{}", self.k - 1);

                },
                CarryType::Max => {
                    tsign = "-";
                    let tmplen = self.graph.node_weight(tid).unwrap().abs_ind.len();

                    tbeg = format!("{}", tmplen);
                    tend = format!("{}$", self.k + tmplen - 1);
                },
            }

            output.push_str( format!("E\t{}\t{}{}\t{}{}\t{}\t{}\t{}\t{}\n",
                                     ei.index(),
                                     sid.index(), ssign,
                                     tid.index(), tsign,
                                     sbeg, send,
                                     tbeg, tend,
                                     ).as_str() );
        } );

        let _ = f.write(&output.into_bytes()[..]);
    }


    fn get_gfa_string(&self) -> String {
        let mut output : String = "VN:Z:2.0\n".to_owned();

        // Given the duality of the graph, to do the exportation of it as GFA, it is enough to assume that e.g. the canonical hashes represent the direct strand.
        // With that, the resulting nodes and edges will represent, by construction, correctly both strands.

        // Add nodes
        self.graph.node_indices().map(|ni| {
            output.push_str( format!("S\t{}\t{}\t*\n", ni.index(), self.k + self.graph.node_weight(ni).unwrap().abs_ind.len() - 1).as_str() );
        } );

        // Add edges
        self.graph.edge_indices().map(|ei| {
            let (sid, tid) = self.graph.edge_endpoints(ei).unwrap();
            let (st, tt)   = self.graph.edge_weight(ei).unwrap().t.get_from_and_to();

            let ssign : &str;
            let tsign : &str;
            let sbeg  : String;
            let send  : String;
            let tbeg  : String;
            let tend  : String;

            match st {
                CarryType::Min => {
                    ssign = "+";

                    let tmplen = self.graph.node_weight(sid).unwrap().abs_ind.len();
                    sbeg = format!("{}",  tmplen);
                    send = format!("{}$", self.k + tmplen - 1);
                },
                CarryType::Max => {
                    ssign = "-";

                    sbeg = format!("{}", 0);
                    send = format!("{}", self.k - 1);
                },
            }

            match tt {
                CarryType::Min => {
                    tsign = "+";
                    tbeg = format!("{}", 0);
                    tend = format!("{}", self.k - 1);

                },
                CarryType::Max => {
                    tsign = "-";
                    let tmplen = self.graph.node_weight(tid).unwrap().abs_ind.len();

                    tbeg = format!("{}", tmplen);
                    tend = format!("{}$", self.k + tmplen - 1);
                },
            }

            output.push_str( format!("E\t{}\t{}{}\t{}{}\t{}\t{}\t{}\t{}\n",
                                     ei.index(),
                                     sid.index(), ssign,
                                     tid.index(), tsign,
                                     sbeg, send,
                                     tbeg, tend,
                                     ).as_str() );
        } );

        output
    }
}



impl Init for PtGraph {
    fn init(edge_count: Option<usize>, node_count: Option<usize>, kl : usize) -> PtGraph {
        let nodes = match node_count {
            Some(n) => n,
            None => 0,
        };
        let edges = match edge_count {
            Some(n) => n,
            None => 0,
        };
        PtGraph {graph : petgraph::stable_graph::StableGraph::with_capacity(nodes, edges),
                 k     : kl,}
    }
}

impl Build for PtGraph {
    fn create_from_map<T: Sized + Init + Build>(k : usize,
                                                    indict : &HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>) -> Self {

        let mut collection = PtGraph::default();
        collection.k = k;
        let mut tmpdict : HashMap::<u64, NodeIndex, BuildHasherDefault<NoHashHasher<u64>>> = HashMap::with_capacity_and_hasher(indict.len(), BuildHasherDefault::default());

        for (h, hi) in indict {
            // First, check if the node exists
            if !tmpdict.contains_key(h) {
                // The node exists, so this must have been added because of some edge. We just have to add the
                // missing edges.
                let tmpid = collection.graph.add_node(NodeStruct {
                                                                counts : hi.borrow().counts,
                                                                abs_ind: vec![*h],
                                                                innerdir: None,
                });
                tmpdict.insert(*h, tmpid);
            }

            // ...and also all of its preceding edges...
            for hpre in hi.borrow().pre.iter() {
                if !tmpdict.contains_key(&hpre.0) {
                    // First, we add the node
                    let tmpid2 = collection.graph.add_node(NodeStruct {
                                                            counts : indict.get(&hpre.0).unwrap().borrow().counts,
                                                            abs_ind: vec![hpre.0],
                                                            innerdir: None,
                    });

                    tmpdict.insert(hpre.0, tmpid2);
                }

                // Then, we add the edge
                if !collection.graph.edges_connecting(*tmpdict.get(&hpre.0).unwrap(), *tmpdict.get(&h).unwrap()).any(|e| e.weight().t == hpre.1) {
                    collection.graph.add_edge(*tmpdict.get(&hpre.0).unwrap(), *tmpdict.get(&h).unwrap(), EmptyEdge{t:hpre.1});
                }

            }
            // ...and the forward ones
            for hpost in hi.borrow().post.iter() {
                if !tmpdict.contains_key(&hpost.0) {
                    // First, we add the node
                    let tmpid2 = collection.graph.add_node(NodeStruct {
                                                            counts : indict.get(&hpost.0).unwrap().borrow().counts,
                                                            abs_ind: vec![hpost.0],
                                                            innerdir: None,
                    });

                    tmpdict.insert(hpost.0, tmpid2);
                }

                // Then, we add the edge
                if !collection.graph.edges_connecting(*tmpdict.get(&h).unwrap(), *tmpdict.get(&hpost.0).unwrap()).any(|e| e.weight().t == hpost.1) {
                    collection.graph.add_edge(*tmpdict.get(&h).unwrap(), *tmpdict.get(&hpost.0).unwrap(), EmptyEdge{t:hpost.1});
                }

            }
        }

        println!("Number of nodes: {}, number of edges: {}", collection.graph.node_count(), collection.graph.edge_count());
        collection
    }
}
