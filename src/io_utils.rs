//! Common helper functions for parsing file input, loading, and setting output
//!
//! The functions are used by a few different subcommands to set correct
//! args to build structs, given the command line input

// use std::error::Error;
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// use regex::Regex;

// use super::QualOpts;
use crate::preprocessing::InputFastx;
// use crate::bit_encoding::UInt;

// use crate::cli::{
//     DEFAULT_KMER, DEFAULT_MINCOUNT, DEFAULT_MINQUAL, DEFAULT_STRAND,
// };



/// Set a buffered stream to write to.
///
/// Either a file (if [`Some`]) or stdout otherwise (if [`None`]).
pub fn set_ostream(oprefix: &Option<String>) -> BufWriter<Box<dyn Write>> {
    let out_writer = match oprefix {
        Some(prefix) => {
            let path = Path::new(prefix);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };
    BufWriter::new(out_writer)
}

/// Obtain a list of input files and names from command line input.
///
/// If `file_list` is provided, read each line as `name\tseq1\tseq2`, where
/// `seq2` is optional, and if present the reverse fastqs. Otherwise, treat
/// as fasta.
///
/// If `seq_files` are provided use [`read_input_fastas`].
pub fn get_input_list(
    file_list: &Option<String>,
    seq_files: &Option<Vec<String>>,
) -> Vec<InputFastx> {
    // Read input
    match file_list {
        Some(files) => {
            let mut input_files: Vec<InputFastx> = Vec::new();
            let f = File::open(files).expect("Unable to open file_list");
            let f = BufReader::new(f);
            for line in f.lines() {
                let line = line.expect("Unable to read line in file_list");
                let fields: Vec<&str> = line.split_whitespace().collect();
                // Should be 2 entries for fasta, 3 for fastq
                let second_file = match fields.len() {
                    0..=1 => {
                        panic!("Unable to parse line in file_list")
                    }
//                     2 => None,
                    2 => panic!("Single-stranded reads are not supported at this moment."),
                    3 => Some(fields[2].to_string()),
                    _ => {
                        panic!("Unable to parse line in file_list")
                    }
                };
                input_files.push((fields[0].to_string(), fields[1].to_string(), second_file));
            }
            input_files
        }
        None => {
            let mut input_files: Vec<InputFastx> = Vec::new();
            let tmpvec = seq_files.as_ref().unwrap();
            input_files.push(("reads".to_owned(), tmpvec[0].clone(), Some(tmpvec[1].clone()) ));
            input_files
        },
        // None => panic!("Single-stranded reads are not supported at this moment."),
    }
}

/// Checks if any input files are fastq
pub fn any_fastq(files: &[InputFastx]) -> bool {
    let mut fastq = false;
    for file in files {
        if file.2.is_some() {
            fastq = true;
            break;
        }
    }
    fastq
}
