//! Graph builder.
use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault, cell::*};

use crate::HashInfoSimple;

/// Custom init function for collections
pub trait Init: Default {
    /// Initialize collection. Arguments are estimated maximum counts of nodes and
    /// edges, as well as type of the input file.
    fn init(_edges_count: Option<usize>, _nodes_count: Option<usize>, k : usize) -> Self {
        Self::default()
    }
}

/// Description of how graph should be built.
pub trait Build: Init {
    /// Create a graph from a map/dictionary structure
    fn create_from_map<T: Sized + Init + Build>(k      : usize,
                                                indict : &HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>) -> Self;
}


