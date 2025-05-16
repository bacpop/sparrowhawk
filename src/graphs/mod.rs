use super::algorithms::builder::Build;
use super::algorithms::collapser::Collapsable;
use super::algorithms::pruner::Prunable;
use pt_graph::{CarryType, NodeIndex};
use std::io::Write;
use pt_graph::EdgeType;

/// Petgraph-based graph for doing DNA de Bruijn's graphs
pub mod pt_graph;

/// Graph's interface.
pub trait Graph
    : Build + Prunable + Collapsable {
    /// Node identifier.
    type NodeIdentifier;
    /// Collection storing nodes which are ambiguous nodes.
    ///
    /// Node is considered ambiguous if this condition holds:
    /// ```(in_degree > 1 || out_degree > 1) || (in_degree == 0 && out_degree >= 1)```
    /// where `in_degree` and `out_degree` are counts of incoming and outgoing
    /// edges.
    type AmbiguousNodes;
    /// Finds all ambiguous in the `Graph`.
    fn get_ambiguous_nodes(&self) -> Self::AmbiguousNodes;
    /// Gets number of outgoing edges for the given node.
    fn get_ambiguous_nodes_bi(&self) -> Self::AmbiguousNodes;
    /// Gets number of outgoing edges for the given node.
    fn node_has_self_loops(&self, node: Self::NodeIdentifier) -> bool;
    /// Gets number of outgoing edges for the given node.
    fn get_good_connections_degree(&self, node: Self::NodeIdentifier) -> usize;
    /// Gets number of outgoing edges for the given node.
    fn get_good_neighbours_bi(&self, node: Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of outgoing edges for the given node.
    fn out_degree(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of outgoing edges for a given node taking into account the canonicality from which you start.
    fn out_degree_bi(&self, nodeid : Self::NodeIdentifier, ty : CarryType) -> usize;
    /// Gets number of outgoing edges for a given node starting from the canonical hash of the node
    fn out_degree_min(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of outgoing edges for a given node starting from the non-canonical hash of the node
    fn out_degree_max(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of incoming edges for a given node.
    fn in_degree(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of incoming edges for a given node taking into account the canonicality from which you start.
    fn in_degree_bi(&self, nodeid : Self::NodeIdentifier, ty : CarryType) -> usize;
    /// Gets number of incoming edges for a given node ending in the canonical hash of the node.
    fn in_degree_min(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of incoming edges for a given node ending in the non-canonical hash of the node.
    fn in_degree_max(&self, nodeid : Self::NodeIdentifier) -> usize;
    /// Gets number of neighbours linked with the provided node that have incoming edges to the provided canonicality.
    fn in_neighbours_bi(&self, node: Self::NodeIdentifier, ty : CarryType) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of neighbours linked with the provided node that have outgoing edges from the provided canonicality.
    fn out_neighbours_bi(&self, nodeid : Self::NodeIdentifier, ty : CarryType) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of neighbours linked with the provided node that have outgoing edges from its canonical hash.
    fn out_neighbours_min(&self, nodeid : Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of neighbours linked with the provided node that have outgoing edges from its non-canonical hash.
    fn out_neighbours_max(&self, nodeid : Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of neighbours linked with the provided node that have incoming edges from its canonical hash.
    fn in_neighbours_min(&self, nodeid : Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)>;
    /// Gets number of neighbours linked with the provided node that have incoming edges from its non-canonical hash.
    fn in_neighbours_max(&self, nodeid : Self::NodeIdentifier) -> Vec<(NodeIndex, EdgeType)>;
    /// Write everything to a graph in DOT format.
    fn write_to_dot<W: Write>(&self, path : &mut W);
    /// Get a String with the graph just pre-collapse in DOT format.
    fn get_dot_string(&self) -> String;
    /// Get the "externals" of a bi-directed DNA de Bruijn graph, i.e. those nodes whose edges only connect them to one only neighbour.
    fn externals_bi(&self) -> Vec<NodeIndex>;
    /// Get the "externals" of a bi-directed DNA de Bruijn graph whose link with their only neighbour is an outgoing edge from the canonical hash of the external.
    fn externals_bi_min(&self, thedir : petgraph::EdgeDirection) -> Vec<NodeIndex>;
    /// Get the "externals" of a bi-directed DNA de Bruijn graph whose link with their only neighbour is an outgoing edge from the non-canonical hash of the external.
    fn externals_bi_max(&self, thedir : petgraph::EdgeDirection) -> Vec<NodeIndex>;
}
