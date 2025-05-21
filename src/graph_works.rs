use std::cell::*;
use std::time::Instant;
use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault};

use rayon::prelude::*;
use std::io::Write;
use super::io_utils::*;

#[cfg(not(feature = "wasm"))]
use needletail::parser::write_fasta;


// use std::process::exit;
extern crate petgraph;
use super::HashInfoSimple;

use crate::graphs::Graph;
use crate::algorithms::collapser::SerializedContigs;
use crate::nthash;
use crate::graphs::pt_graph::EdgeType;

use crate::bit_encoding::rc_base;


/// Get backwards neighbours, i.e. incoming edges to either the canonical or non-canonical hashes
pub fn check_bkg( // backwards here mean INCOMING edges, whether from the rev.-comp. strand, or the direct one
    hc:         u64,
    hnc:        u64,
    k:          usize,
    bases:      u8,
    thedict:    &HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>,
    maxmindict: &HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>,
) -> Vec<(u64, EdgeType)>
{
    let mut outvec = Vec::new();
    // Let's start first with the canonical one
    let thecbase  = bases & 3;
    let thencbase = (bases >> 2) & 3;
    for i in 0..4 {             //ACTG, in that order
        let tmphashc = (nthash::swapbits033(hc ^ nthash::HASH_LOOKUP[thecbase as usize]
                                               ^ (nthash::MS_TAB_31L[(i as usize * 31) + (k % 31)]
                                                | nthash::MS_TAB_33R[(i as usize) * 33 + (k % 33)])
                        )).rotate_right(1u32);

        if thedict.contains_key(&tmphashc) {
            outvec.push((tmphashc, EdgeType::MinToMin));
        } else {
            let poth = maxmindict.get(&tmphashc);
            if poth.is_some_and(|x| thedict.contains_key(x)) {
                outvec.push((*poth.unwrap(), EdgeType::MaxToMin));
            }
        }

        let mut tmphashnc = hnc
            ^ (nthash::MS_TAB_31L[(rc_base(i) as usize  * 31) + (k % 31)]
            |  nthash::MS_TAB_33R[(rc_base(i) as usize) * 33  + (k % 33)]);
        tmphashnc ^= nthash::RC_HASH_LOOKUP[thencbase as usize];
        tmphashnc = tmphashnc.rotate_right(1_u32);
        tmphashnc = nthash::swapbits3263(tmphashnc);


        if thedict.contains_key(&tmphashnc) {
            outvec.push((tmphashnc, EdgeType::MinToMax));
        } else {
            let poth = maxmindict.get(&tmphashnc);
            if poth.is_some_and(|x| thedict.contains_key(x)) {
                outvec.push((*poth.unwrap(), EdgeType::MaxToMax));
            }
        }
    }

    outvec
}


/// Get forward neighbours, i.e. outgoing edges from either the canonical or non-canonical hashes
pub fn check_fwd( // Here FORWARD means OUTGOING
    hc:      u64,
    hnc:     u64,
    k:       usize,
    bases:   u8,
    thedict: &HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>,
    maxmindict: &HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>,
) -> Vec<(u64, EdgeType)>
{
    let mut outvec = Vec::new();

    // Let's start first with the canonical one
    let thecbase  = (bases >> 2) & 3;
    let thencbase  = bases & 3;
    for i in 0..4 {             //ACTG, in that order
        let mut tmphashc = hc.rotate_left(1);
        tmphashc =  nthash::swapbits033(tmphashc);
        tmphashc ^= nthash::HASH_LOOKUP[i as usize];
        tmphashc ^= nthash::MS_TAB_31L[(thecbase as usize * 31) + (k % 31)]
                  | nthash::MS_TAB_33R[(thecbase as usize) * 33 + (k % 33)];


        if thedict.contains_key(&tmphashc) {
            outvec.push((tmphashc, EdgeType::MinToMin));
        } else {
            let poth = maxmindict.get(&tmphashc);
            if poth.is_some_and(|x| thedict.contains_key(x)) {
                outvec.push((*poth.unwrap(), EdgeType::MinToMax));
            }
        }

        let tmphashnc = (nthash::swapbits3263(hnc)).rotate_left(1u32)
            ^ nthash::RC_HASH_LOOKUP[i as usize]
            ^ (nthash::MS_TAB_31L[(rc_base(thencbase) as usize * 31) + (k % 31)]
             | nthash::MS_TAB_33R[(rc_base(thencbase) as usize) * 33 + (k % 33)]);


        if thedict.contains_key(&tmphashnc) {
            outvec.push((tmphashnc, EdgeType::MaxToMin));
        } else {
            let poth = maxmindict.get(&tmphashnc);
            if poth.is_some_and(|x| thedict.contains_key(x)) {
                outvec.push((*poth.unwrap(), EdgeType::MaxToMax));
            }
        }
    }

    outvec
}


////////////////////////////////////////////////////////////////////////
/// Output from the assembler.
#[derive(Default)]
pub struct Contigs {
    /// Serialized contigs.
    pub serialized_contigs: SerializedContigs,

    /// Sequence of the contigs.
    pub contig_sequences: Option<Vec<Vec<u8>>>,
}


impl Contigs {
    /// Create new `Contigs`.
    pub fn new(serialized: SerializedContigs) -> Contigs {
        Contigs {
            serialized_contigs: serialized,
            contig_sequences: None,
        }
    }


    /// Temporal and historical function to simplify contigs. To be removed in the future
    pub fn shrink(&mut self) {
        log::warn!("Contigs are going to be shrunk!");

        for ic in 0..self.serialized_contigs.len() {
            let mut tmpv = self.serialized_contigs[ic][0].abs_ind.clone();
            let contiglen = self.serialized_contigs[ic].len();
            for j in 1..contiglen {
                tmpv.extend(self.serialized_contigs[ic][j].abs_ind.clone());
            }
            self.serialized_contigs[ic].first_mut().unwrap().abs_ind = tmpv;
            self.serialized_contigs[ic].drain(1..);
        }
    }


    /// Save contigs in a file
    #[cfg(not(feature = "wasm"))]
    pub fn write_fasta<W: Write>(&self, f: &mut W) {
//         self.print_coverage_stats();
        for i in 0..self.contig_sequences.as_ref().unwrap().len() {
            let _ = write_fasta(
                i.to_string().as_bytes(),
                &self.contig_sequences.as_ref().unwrap()[i][..],
                f,
                needletail::parser::LineEnding::Unix,
            );
        }
    }
}


/// Public API for assemblers.
pub trait Assemble {
    /// Assembles given data using specified `Graph` and writes results into the output file.
    fn assemble<G: Graph>(k : usize, indict : &mut HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>,
                          maxminsize : &mut HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, timevec : &mut Vec<Instant>,
                          path : Option<&String>) -> (Contigs, String, String, String);
}


///////////////////////////////////////////////////////////
/// Basic assembler.
pub struct BasicAsm {}


impl Assemble for BasicAsm {
    fn assemble<G: Graph>(k    : usize,
                        indict     : &mut HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>,
                        maxmindict : &mut HashMap::<u64, u64,                     BuildHasherDefault<NoHashHasher<u64>>>,
                        timevec    : &mut Vec<Instant>, path : Option<&String>) -> (Contigs, String, String, String) {
        log::info!("Starting assembler!");

        // FIRST: iterate over all k-mers, check the existance of forwards/backwards neighbours in the dictionary.
        let mut i = 0;
        let mut ialone = 0;
        let mut nedges = 0;
        indict.iter().for_each(|(h, hi)| {
            let mut himutref = hi.borrow_mut();

            himutref.pre  = check_bkg(*h, himutref.hnc, k, himutref.b, indict, maxmindict);
            himutref.post = check_fwd(*h, himutref.hnc, k, himutref.b, indict, maxmindict);
            let tmpnedges = himutref.pre.len() + himutref.post.len();
            nedges += tmpnedges;
            i += 1;
            if tmpnedges == 0 { ialone += 1};
        });

//         drop(maxmindict);
        println!("Prop. of alone kmers: {:.1} %", (ialone as f64) / (i as f64) * 100.0);
        println!("Number of edges {}", (nedges as f64) / (2 as f64));

        // timevec.push(Instant::now());
        // log::info!("Neighbours searched for in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());

        // indict.iter().for_each(|(h, hi)| {
        //     let himutref = hi.borrow();
        //
        //     // Check first previous neighbours:
        //     for ipre in himutref.pre.iter() {
        //         if !indict.get(&ipre.0).unwrap().borrow().pre.contains(&(*h, ipre.1.rev())) || !indict.get(&ipre.0).unwrap().borrow().post.contains(&(*h, ipre.1)) {
        //             println!("HEY");
        //         }
        //     }
        //     for ipost in himutref.post.iter() {
        //         if !indict.get(&ipost.0).unwrap().borrow().post.contains(&(*h, ipost.1.rev())) || !indict.get(&ipost.0).unwrap().borrow().pre.contains(&(*h, ipost.1)) {
        //             println!("HEY");
        //         }
        //     }
        // });

        let ptgraph = G::create_from_map::<G>(k, indict);

        assemble_with_bi_graph(ptgraph, timevec, path)
    }
}


/// Assemble a bidirected DNA de Bruijn graph
fn assemble_with_bi_graph<G: Graph>(mut ptgraph: G, timevec : &mut Vec<Instant>, path_ : Option<&String>)
-> (Contigs, String, String, String) {

    // log::info!("Saving graph (pre-shrink w/o one-node contigs) as DOT file...");
    // if path_.is_some() {
    //     let mut wbuf = set_ostream(&Some(path_.unwrap().clone().replace(".dot", "_preshrink.dot")));
    //     ptgraph.write_to_dot(&mut wbuf);
    // }
    // log::info!("Done.");

    log::info!("Starting shrinkage and pruning of the graph");

    log::info!("Removing self-loops (temporal restriction)");
    ptgraph.remove_self_loops();

    // let minnts = 100; // independent of this value, the minimum number of nts will be always k
    // let limit = max(0, minnts - k + 1);

    let mut didanyofusdoanything = true;
    while didanyofusdoanything {
        let bools = ptgraph.shrink();
        // let bools = false;
        let boolr = ptgraph.remove_dead_paths();
        // let boolr = false;
//             println!("{} {}", bools, boolr);
        didanyofusdoanything = bools || boolr;
        // break;
    }
    log::info!("Shrinkage and pruning finished");


    log::info!("Saving graph (post-shrink, pre-collapse, w/o one-node contigs) as DOT file...");
    if path_.is_some() {
        let mut wbuf = set_ostream(&Some(path_.unwrap().clone()));
        ptgraph.write_to_dot(&mut wbuf);
    }
    log::info!("Done.");

    let outdot  = ptgraph.get_dot_string();
    let outgfa  = ptgraph.get_gfa_string();
    let outgfa2 = ptgraph.get_gfa2_string();

    let serialized_contigs = ptgraph.collapse(path_);
    log::info!("I created {} contigs", serialized_contigs.len());
    let mut contigs = Contigs::new(serialized_contigs);

    // TEMPORAL RESTRICTION, WIP
    contigs.shrink();

    (contigs, outdot, outgfa, outgfa2)
}
