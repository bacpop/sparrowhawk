//! Efficient genome assembler for small genomes in Rust
#![warn(missing_docs)]
use std::fmt;
use std::time::Instant;

extern crate num_cpus;


/// Construction, assembly, shrinkage, pruning, and collapse of DNA de Bruijn graphs
pub mod graph_works;

/// Preprocessing functions of the reads & k-mers
pub mod preprocessing;

/// Declarations and definitions for encode our k-mers efficiently in memory
pub mod bit_encoding;

/// A helper class to obtain k-mers from reads derived from Ska2
pub mod kmer;

/// An implementation of ntHash, based on ntHash 2
pub mod nthash;

/// Contains functions to store the output of the program
pub mod save_functions;

/// Contains the graph definitions
pub mod graphs;

/// Contains different traits that implement various algorithms
pub mod algorithms;


use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault, cell::*};
use crate::graphs::pt_graph::EdgeType;

use crate::graphs::pt_graph::PtGraph;
use crate::graph_works::Assemble;
use bit_encoding::{U256, U512};

pub mod cli;
use crate::cli::*;

pub mod io_utils;
use crate::io_utils::*;


#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;
#[cfg(feature = "wasm")]
use wasm_bindgen_file_reader::WebSysFile;
#[cfg(feature = "wasm")]
extern crate console_error_panic_hook;
#[cfg(feature = "wasm")]
pub mod fastx_wasm;
#[cfg(feature = "wasm")]
use crate::graph_works::Contigs;
#[cfg(feature = "wasm")]
use json;
#[cfg(feature = "wasm")]
use crate::bit_encoding::UInt;

// use std::process::exit;


#[cfg(feature = "wasm")]
pub fn loG(text : &str, typ : Option<&str>) {
    if typ.is_some() {
        log((String::from("Sparrowhawk::") + typ.unwrap() + "::" + text).as_str());
    } else {
        log(text);
    }
}

#[cfg(not(feature = "wasm"))]
pub fn loG(text : &str, typ : Option<&str>) {
    if let Some(realtyp) = typ {
        if realtyp == "info" {
            log::info!("{}", text);
        } else if realtyp == "debug" {
            log::debug!("{}", text);
        } else if realtyp == "trace" {
            log::trace!("{}", text);
        } else if realtyp == "warn" {
            log::warn!("{}", text);
        } else if realtyp == "error" {
            log::error!("{}", text);
        } else {
            println!("{}", text);
        }
    } else {
        println!("{}", text);
    }
}


/// Index type for both nodes and edges in the graph/gir.
pub type Idx = usize;


/// Type for representing weight of the `Edge`.
pub type EdgeWeight = u16;


/// Struct that contains the basic information for one k-mer
pub struct HashInfoSimple {
    /// maximum hash
    pub hnc: u64,
    /// First and last bases
    pub b: u8,
    /// found neighbours, if any, previous to this kmer
    pub pre: Vec<(u64, EdgeType)>,
    /// found neighbours, if any, posterior to this kmer
    pub post: Vec<(u64, EdgeType)>,
    /// Counts associated to this kmer
    pub counts: u16,
}


/// Quality filtering options for FASTQ files
pub struct QualOpts {
    /// Minimum k-mer count across reads to be added
    pub min_count: u16,
    /// Minimum base quality to be added
    pub min_qual: u8,
}


impl fmt::Display for QualOpts {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "min count: {}; minimum quality {} ({});",
            self.min_count,
            self.min_qual,
            (self.min_qual + 33) as char,
        )
    }
}


#[doc(hidden)]
#[cfg(not(feature = "wasm"))]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Trace).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    log::info!("Starting program!");
    eprintln!("Sparrowhawk");
    let mut timevec = Vec::new();
    timevec.push(Instant::now());
    match &args.command {
        Commands::Build {
            seq_files,
            file_list,
            output,
            k,
            single_strand,
            min_count,
            min_qual,
            threads,
        } => {
            check_threads(*threads);

            // Read input
            // let input_files = get_input_list(file_list, seq_files);
            let input_files = get_input_list(file_list);
            let quality = QualOpts {
                min_count: *min_count,
                min_qual: *min_qual,
            };

            // Build, merge
            // let rc = !*single_strand;

            log::info!("Checking requested threads and creating pool if needed");
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();
            log::info!("Beginning processing");
            timevec.push(Instant::now());


            let mut preprocessed_data : HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>;
            let theseq : Vec<u64>;
            let mut maxmindict : HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>;

            if *k % 2 == 0 {
                panic!("Support for even k-mer lengths not implemented");

            } else if *k < 3 {
                panic!("kmer length too small (min. 3)");

            } else if *k <= 32 {
                log::info!("k={}: using 64-bit representation", *k);
                let thedict : HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>;
                (preprocessed_data, theseq, thedict, maxmindict) = preprocessing::preprocessing_gpulike_with_dict_and_seq::<u64>(&input_files, *k, &quality, &mut timevec);
                drop(theseq);
                let outgraph = output.to_string().clone().replace(".fasta", "") + ".dot";
                let mut contigs = graph_works::BasicAsm::assemble::<PtGraph>(*k, &mut preprocessed_data, &mut maxmindict, &mut timevec, Some(&outgraph));

                // Save as fasta
                save_functions::save_as_fasta::<u64>(&mut contigs, &thedict, *k, output); // FASTA file(s)

            } else if *k <= 64 {
                log::info!("k={}: using 128-bit representation", *k);
                let thedict : HashMap::<u64, u128, BuildHasherDefault<NoHashHasher<u64>>>;
                (preprocessed_data, theseq, thedict, maxmindict) = preprocessing::preprocessing_gpulike_with_dict_and_seq::<u128>(&input_files, *k, &quality, &mut timevec);
                drop(theseq);

                let outgraph = output.to_string().clone().replace(".fasta", "") + ".dot";
                let mut contigs = graph_works::BasicAsm::assemble::<PtGraph>(*k, &mut preprocessed_data, &mut maxmindict, &mut timevec, Some(&outgraph));

                // Save as fasta
                save_functions::save_as_fasta::<u128>(&mut contigs, &thedict, *k, output); // FASTA file(s)

            } else if *k <= 128 {
                log::info!("k={}: using 256-bit representation", *k);

                let thedict : HashMap::<u64, U256, BuildHasherDefault<NoHashHasher<u64>>>;
                (preprocessed_data, theseq, thedict, maxmindict) = preprocessing::preprocessing_gpulike_with_dict_and_seq::<U256>(&input_files, *k, &quality, &mut timevec);
                drop(theseq);

                let outgraph = output.to_string().clone().replace(".fasta", "") + ".dot";
                let mut contigs = graph_works::BasicAsm::assemble::<PtGraph>(*k, &mut preprocessed_data, &mut maxmindict, &mut timevec, Some(&outgraph));

                // Save as fasta
                save_functions::save_as_fasta::<U256>(&mut contigs, &thedict, *k, output); // FASTA file(s)

            } else if *k <= 256 {
                log::info!("k={}: using 512-bit representation", *k);

                let thedict : HashMap::<u64, U512, BuildHasherDefault<NoHashHasher<u64>>>;
                (preprocessed_data, theseq, thedict, maxmindict) = preprocessing::preprocessing_gpulike_with_dict_and_seq::<U512>(&input_files, *k, &quality, &mut timevec);
                drop(theseq);

                let outgraph = output.to_string().clone().replace(".fasta", "") + ".dot";
                let mut contigs = graph_works::BasicAsm::assemble::<PtGraph>(*k, &mut preprocessed_data, &mut maxmindict, &mut timevec, Some(&outgraph));

                // Save as fasta
                save_functions::save_as_fasta::<U512>(&mut contigs, &thedict, *k, output); // FASTA file(s)

            } else {
                panic!("kmer length larger than 256 currently not supported.");
            }
        }
    }

    timevec.push(Instant::now());

    eprintln!("Sparrowhawk done in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());
    log::info!("Finishing program!");
}


#[cfg(feature = "wasm")]
pub fn main() {
    panic!("You've compiled Sparrowhawk for WebAssembly support, you cannot use it as a normal binary anymore!");
}


// ===================================== WebAssembly stuff follows

#[cfg(feature = "wasm")]
#[wasm_bindgen]
extern {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}


#[cfg(feature = "wasm")]
#[wasm_bindgen]
pub struct AssemblyHelper {
    k                 : usize,
    preprocessed_data : HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>,
    maxmindict        : HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>,
    seqdict64         : Option<HashMap::<u64, u64,  BuildHasherDefault<NoHashHasher<u64>>>>,
    seqdict128        : Option<HashMap::<u64, u128, BuildHasherDefault<NoHashHasher<u64>>>>,
    seqdict256        : Option<HashMap::<u64, U256, BuildHasherDefault<NoHashHasher<u64>>>>,
    seqdict512        : Option<HashMap::<u64, U512, BuildHasherDefault<NoHashHasher<u64>>>>,
    histovec          : Vec<u16>,
    contigs           : Contigs,
    outfasta          : String,
}

#[cfg(feature = "wasm")]
#[wasm_bindgen]
impl AssemblyHelper {
    pub fn new(file1 : web_sys::File, file2 : web_sys::File, k : usize, verbose : bool, min_count : u16, min_qual : u8, csize : usize) -> Self {
        if cfg!(debug_assertions) {
            init_panic_hook();
        }

        let mut wf1 = WebSysFile::new(file1);
        let mut wf2 = WebSysFile::new(file2);

        // Read input
        let quality = QualOpts {
            min_count: min_count,
            min_qual:  min_qual,
        };

        // loG("Checking requested threads and creating pool if needed", Some("info"));
        // rayon::ThreadPoolBuilder::new()
        //     .num_threads(8)
        //     .build_global()
        //     .unwrap();
        loG("Beginning processing", Some("info"));

        let mut preprocessed_data : HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>;
        let mut maxmindict : HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>;
        let histovalues : Vec<u16>;
        let mut thedict64  = None;
        let mut thedict128 = None;
        let mut thedict256 = None;
        let mut thedict512 = None;

        if k % 2 == 0 {
            panic!("Support for even k-mer lengths not implemented");

        } else if k < 3 {
            panic!("kmer length too small (min. 3)");

        } else if k <= 32 {
            loG(format!("k={}: using 64-bit representation", k).as_str(), Some("info"));

            (preprocessed_data, thedict64, maxmindict, histovalues) = preprocessing::preprocessing_for_wasm::<u64>(&mut wf1, &mut wf2, k, &quality, csize);

            loG("Preprocessing done!", Some("info"));

        } else if k <= 64 {
            loG(format!("k={}: using 128-bit representation", k).as_str(), Some("info"));

            (preprocessed_data, thedict128, maxmindict, histovalues) = preprocessing::preprocessing_for_wasm::<u128>(&mut wf1, &mut wf2, k, &quality, csize);

            loG("Preprocessing done!", Some("info"));

        } else if k <= 128 {
            loG(format!("k={}: using 256-bit representation", k).as_str(), Some("info"));

            (preprocessed_data, thedict256, maxmindict, histovalues) = preprocessing::preprocessing_for_wasm::<U256>(&mut wf1, &mut wf2, k, &quality, csize);

            loG("Preprocessing done!", Some("info"));
        } else if k <= 256 {
            loG(format!("k={}: using 512-bit representation", k).as_str(), Some("info"));

            (preprocessed_data, thedict512, maxmindict, histovalues) = preprocessing::preprocessing_for_wasm::<U512>(&mut wf1, &mut wf2, k, &quality, csize);

            loG("Preprocessing done!", Some("info"));
        } else {
            panic!("kmer length larger than 256 currently not supported.");
        }

        // loG("Sparrowhawk done!", Some("info"));

        Self {
            k : k,
            preprocessed_data : preprocessed_data,
            maxmindict        : maxmindict,
            seqdict64         : thedict64,
            seqdict128        : thedict128,
            seqdict256        : thedict256,
            seqdict512        : thedict512,
            histovec          : histovalues,
            contigs           : Contigs::default(),
            outfasta          : "".to_owned(),
        }
    }


    pub fn assemble(&mut self) {
        loG("Starting assembly...", Some("info"));
        let mut outcontigs = graph_works::BasicAsm::assemble::<PtGraph>(self.k, &mut self.preprocessed_data, &mut self.maxmindict, &mut Vec::new(), None);

        loG("Assembly done!", Some("info"));

        let mut outfasta   : String;

        if self.k % 2 == 0 {
            panic!("Support for even k-mer lengths not implemented");

        } else if self.k < 3 {
            panic!("kmer length too small (min. 3)");

        } else if self.k <= 32 {
            outfasta = save_functions::save_as_fasta_for_wasm::<u64>(&mut outcontigs, self.seqdict64.as_ref().unwrap(), self.k);
        } else if self.k <= 64 {
            outfasta = save_functions::save_as_fasta_for_wasm::<u128>(&mut outcontigs, self.seqdict128.as_ref().unwrap(), self.k);
        } else if self.k <= 128 {
            outfasta = save_functions::save_as_fasta_for_wasm::<U256>(&mut outcontigs, self.seqdict256.as_ref().unwrap(), self.k);
        } else if self.k <= 256 {
            outfasta = save_functions::save_as_fasta_for_wasm::<U512>(&mut outcontigs, self.seqdict512.as_ref().unwrap(), self.k);
        } else {
            panic!("kmer length larger than 256 currently not supported.");
        }

        loG("Sparrowhawk done!", Some("info"));

        self.contigs  = outcontigs;
        self.outfasta = outfasta;
    }


    pub fn get_assembly(&self) -> String {
        let mut results = json::JsonValue::new_array();

        results["outfasta"] = json::JsonValue::String(self.outfasta.clone());
        results["ncontigs"] = json::JsonValue::Number(self.contigs.contig_sequences.as_ref().unwrap().len().into());

        return results.dump();
    }


    pub fn get_preprocessing_info(&self) -> String {
        let mut results = json::JsonValue::new_array();

        loG(format!("{} {}", self.preprocessed_data.len(), self.histovec.len()).as_str(), Some("info"));

        results["nkmers"] = json::JsonValue::Number(self.preprocessed_data.len().into());
        results["histo"]  = json::JsonValue::Array(self.histovec.iter().map(|x| json::JsonValue::Number((*x).into())).collect());

        loG(results.dump().as_str(), Some("debug"));

        return results.dump();
    }
}
