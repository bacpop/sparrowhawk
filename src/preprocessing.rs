//! Some docs should be here

use core::panic;

#[cfg(not(feature = "wasm"))]
use std::{
    time::Instant,
    path::PathBuf,
};

use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault, cell::*, cmp::Ordering};

use rayon::prelude::*;

#[cfg(not(feature = "wasm"))]
use needletail::parse_fastx_file;

use std::process::exit;

#[cfg(not(feature = "wasm"))]
use plotters::prelude::*;


use super::QualOpts;
use super::HashInfoSimple;

#[cfg(not(feature = "wasm"))]
use super::bit_encoding::{encode_base, rc_base};

#[cfg(not(feature = "wasm"))]
use crate::io_utils::any_fastq;

use crate::bit_encoding::UInt;
use crate::kmer::Kmer;
use crate::logw;
use crate::spectrum_fitter::SpectrumFitter;
use crate::bloom_filter::KmerFilter;


/// Tuple for name and fasta or paired fastq input
pub type InputFastx = (String, String, Option<String>);

// #[cfg(feature = "wasm")]
// use wasm_bindgen::prelude::*;
#[cfg(feature = "wasm")]
use wasm_bindgen_file_reader::WebSysFile;
#[cfg(feature = "wasm")]
use crate::fastx_wasm::open_fastq;
#[cfg(feature = "wasm")]
use seq_io::fastq::Record;

// For the fitting, we'll use actually MAXSIZEHISTO - 1
const MAXSIZEHISTO : usize = 500;

// =====================================================================================================


#[cfg(not(feature = "wasm"))]
fn put_these_nts_into_an_efficient_vector(charseq : &[u8], compseq : &mut Vec<u64>, occ : u8) {
    let mut tmpu64 : u64 = 0;
    let mut tmpind : u8  = 0;

    if occ != 0 {
        tmpind = occ;
        tmpu64 = compseq.pop().unwrap();
    }
//     log::debug!("{}", tmpind);

    for nt in charseq {
//         log::debug!("\n{:#010b}\n{}\n{:#066b}\n{:#066b}", *nt, tmpind, tmpu64, (encode_base(*nt) as u64));
        tmpu64 <<= 2;
        tmpu64 |= encode_base(*nt) as u64;
//         log::debug!("{:#066b}", tmpu64);
        if tmpind == 31 {
            compseq.push(tmpu64);
            tmpu64 = 0;
            tmpind = 0;
        } else {
            tmpind += 1;
        }
//         log::debug!("{}", tmpind);
    }

    if tmpind != 0 {
        compseq.push(tmpu64);
    }
}

#[cfg(not(feature = "wasm"))]
fn put_these_nts_into_an_efficient_vector_rc(charseq : &[u8], compseq : &mut Vec<u64>, occ : u8) {
    let mut tmpu64 : u64 = 0;
    let mut tmpind : u8  = 0;

    if occ != 0 {
        tmpind = occ;
        tmpu64 = compseq.pop().unwrap();
    }
//     log::debug!("{}", tmpind);
    for nt in charseq.iter().rev() {
//         log::debug!("\n{:#010b}\n{}\n{:#066b}\n{:#066b}", *nt, tmpind, tmpu64, (rc_base(encode_base(*nt)) as u64));
        tmpu64 <<= 2;
        tmpu64 |= rc_base(encode_base(*nt)) as u64;
//         log::debug!("{:#066b}", tmpu64);
        if tmpind == 31 {
            compseq.push(tmpu64);
            tmpu64 = 0;
            tmpind = 0;
        } else {
            tmpind += 1;
        }
    }

    if tmpind != 0 {
        compseq.push(tmpu64);
    }
}



#[cfg(not(feature = "wasm"))]
fn get_kmers_from_both_files_and_the_dict_and_the_seq<IntT>(
    filename1: &str,
    filename2: &str,                    // Will get the "rc" automatically enabled
    k:         usize,
    qual:      &QualOpts,
    outvec:    &mut Vec<(u64, u64, u8)>,
) -> (Vec<u64>, HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    log::info!("Getting kmers from file {filename1}. Creating reader...");
    let mut reader =
        parse_fastx_file(filename1).unwrap_or_else(|_| panic!("Invalid path/file: {filename1}"));

    log::info!("Entering while loop...");

//     let maxkmers = 200;
//     let mut numkmers = 0;
    let mut ncols : usize = 0;
    let mut itrecord : u32 = 0;                                 // We're using it to add the previous indexes!
    let mut theseq   : Vec<u64> = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        put_these_nts_into_an_efficient_vector(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);

        let rl = seqrec.num_bases();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
            // &itrecord,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            // outdict.entry(hc).or_insert(km);
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                log::debug!("\n\t- COLLISIONS 1 !!! Hash: {:?}", hc);
                log::debug!("{:#0258b}\n{:#0258b}", *testkm, km);
                ncols += 1;
            }
            // } else {
            //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                // outdict.entry(hc).or_insert(km);
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    log::debug!("\n\t- COLLISIONS 2 !!! Hash: {:?}", hc);
                    log::debug!("{:#0258b}\n{:#0258b}", *testkm, km);
                    ncols += 1;
                }
                // } else {
                //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
                // }
                minmaxdict.entry(hnc).or_insert(hc);

//                 numkmers += 1;
//                 if numkmers >= maxkmers {break};
            }
        }
//         if numkmers >= maxkmers {itrecord += rl as u32;break};
        itrecord += rl as u32;
    }
    log::info!("Finished getting kmers from file {filename1}. Starting with {filename2}");
//     numkmers = 0;

    let mut reader =
        parse_fastx_file(filename2).unwrap_or_else(|_| panic!("Invalid path/file: {filename2}"));

//     log::debug!("MITAD Length of seq. vec.: {}, total length of first files: {}, THING {}, number of recs: {}", theseq.len(), itrecord, (itrecord % 32) as u8, realit);
    // Memory usage optimisations
    let cseq = theseq.capacity();
    let lseq = theseq.len();
    if cseq < 2 * lseq {
        theseq.reserve_exact(2 * lseq);
    } else {
        theseq.shrink_to(2 * lseq);
    }
    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
//         log::debug!("{}", seqrec.num_bases());
        let rl = seqrec.num_bases();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
            // &itrecord,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            // outdict.entry(hc).or_insert(km);
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                log::debug!("\n\t- COLLISIONS 3 !!! Hash: {:?}", hc);
                log::debug!("{:#0258b}\n{:#0258b}", *testkm, km);
                ncols += 1;
            }
            // } else {
            //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                // outdict.entry(hc).or_insert(km);
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    log::debug!("\n\t- COLLISIONS 4 !!! Hash: {:?}", hc);
                    log::debug!("{:#0258b}\n{:#0258b}", *testkm, km);
                    ncols += 1;
                }
                // } else {
                //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
                // }
                minmaxdict.entry(hnc).or_insert(hc);

//                 numkmers += 1;
//                 if numkmers >= maxkmers {break};
            }
        }
//         if numkmers >= maxkmers {itrecord += rl as u32;break};
        itrecord += rl as u32;
    }

    if (itrecord % 32) != 0 {
        let mut tmpu64 = theseq.pop().unwrap();
        tmpu64 <<= 2 * (32 - itrecord % 32);
//         log::debug!("{:#066b}", tmpu64);
        theseq.push(tmpu64);
    }

    log::info!("Finished getting kmers from file {filename2}");
    log::debug!("Length of seq. vec.: {}, total length of both files: {}", theseq.len(), itrecord);
    log::debug!("k | Number of collisions =+=+ {} {}", k, ncols);


    (theseq, outdict, minmaxdict)
}



#[cfg(feature = "wasm")]
fn get_kmers_from_both_files_wasm<IntT>(
    file1:    &mut WebSysFile,
    file2:    &mut WebSysFile,
    k:         usize,
    qual:      &QualOpts,
    outvec:    &mut Vec<(u64, u64, u8)>,
) -> (HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    logw("Getting kmers from first file. Creating reader...", Some("info"));
    let mut reader = open_fastq(file1);

    logw("Entering while loop...", Some("info"));

//     let maxkmers = 200;
//     let mut numkmers = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            // let testkm = outdict.entry(hc).or_insert(km);
            // if *testkm != km {
            //     logw(format!("\n\t- COLLISIONS 1 !!! Hash: {:?}", hc).as_str(), Some("warn"));
            //     logw(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
            // }
            // } else {
            //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                // let testkm = outdict.entry(hc).or_insert(km);
                // if *testkm != km {
                //     logw(format!("\n\t- COLLISIONS 2 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                //     logw(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
                // }
                // } else {
                //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
                // }
                minmaxdict.entry(hnc).or_insert(hc);

//                 numkmers += 1;
//                 if numkmers >= maxkmers {break};
            }
        }
    }
    logw("Finished getting kmers from first file. Starting with the second...", Some("info"));
//     numkmers = 0;

    let mut reader = open_fastq(file2);

//     log::debug!("MITAD Length of seq. vec.: {}, total length of first files: {}, THING {}, number of recs: {}", theseq.len(), itrecord, (itrecord % 32) as u8, realit);
    // Memory usage optimisations
    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            // let testkm = outdict.entry(hc).or_insert(km);
            // if *testkm != km {
            //     logw(format!("\n\t- COLLISIONS 3 !!! Hash: {:?}", hc).as_str(), Some("warn"));
            //     logw(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
            // }
            // } else {
            //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                // let testkm = outdict.entry(hc).or_insert(km);
                // if *testkm != km {
                //     logw(format!("\n\t- COLLISIONS 4 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                //     logw(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
                // }
                // } else {
                //     log::debug!("\n\t\t- NOT COLLISIONS!!!");
                // }
                minmaxdict.entry(hnc).or_insert(hc);

//                 numkmers += 1;
//                 if numkmers >= maxkmers {break};
            }
        }
//         if numkmers >= maxkmers {itrecord += rl as u32;break};
        // itrecord += rl as u32;
    }

    logw("Finished getting kmers from the second file", Some("info"));


    (outdict, minmaxdict)
}


#[cfg(feature = "wasm")]
fn chunked_processing_wasm<IntT>(
    file1:    &mut WebSysFile,
    file2:    &mut WebSysFile,
    k:         usize,
    qual:      &QualOpts,
    outvec:    &mut Vec<(u64, u64, u8)>,
    csize:     usize,
    do_fit:    bool,
) -> (HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u32>, u16)
where
    IntT: for<'a> UInt<'a>,
{
    logw("Getting kmers from first file. Creating reader...", Some("info"));

    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    let mut themap     = HashMap::with_hasher(BuildHasherDefault::default());
    let mut countmap : HashMap::<u64, (u16, u64, u8), BuildHasherDefault<NoHashHasher<u64>>> = HashMap::with_hasher(BuildHasherDefault::default());

    let mut reader = open_fastq(file1);
    let mut histovec : Vec<u32> = vec![0; MAXSIZEHISTO];

    logw("Entering while loop...", Some("info"));

    let mut i_record = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
        }

        i_record += 1;
        if i_record >= csize {
            // Processssssss! And reset.
            if !outvec.is_empty() {
                logw("Processing chunk. Sorting k-mers...", Some("info"));
                outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                logw("k-mers sorted. Counting k-mers...", Some("info"));
                // Then, do a counting of everything and save the results in a dictionary and return it

                update_countmap(&outvec, &mut countmap);
            }

            // Reset
            outvec.clear();
            i_record = 0;
        }

    }
    logw("Finished getting kmers from first file. Starting with the second...", Some("info"));

    let mut reader = open_fastq(file2);

    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        // put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
        }

        i_record += 1;
        if i_record >= csize {
            // Processssssss! And reset.
            if !outvec.is_empty() {
                logw("Processing chunk. Sorting k-mers...", Some("info"));
                outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                logw("k-mers sorted. Counting k-mers...", Some("info"));
                // Then, do a counting of everything and save the results in a dictionary and return it

                update_countmap(&outvec, &mut countmap);
            }

            // Reset
            outvec.clear();
            i_record = 0;
        }
    }

    if i_record > 0 {
        // Processssssss! And reset.
        if !outvec.is_empty() {
            logw("Processing last chunk. Sorting k-mers...", Some("info"));
            outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            logw("k-mers sorted. Counting k-mers...", Some("info"));
            // Then, do a counting of everything and save the results in a dictionary and return it

            update_countmap(&outvec, &mut countmap);
        }
        // Reset
        outvec.clear();
    }
    logw("Finished getting kmers from the second file", Some("info"));
    logw("Filtering...", Some("info"));

    // Now, get themap, histovec, and filter outdict and minmaxdict
    countmap.shrink_to_fit();
    let mut minc : u16 = qual.min_count;


    // This can be optimised. also better written: I had to repeat the code for the retains, to try to improve slightly the running time in
    // case no autofitting is requested. In any case, it could be improved in the future.
    if do_fit {
        for (_, tup) in countmap.iter() {
            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
        }

        // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
        // large (and so that we can detect it).
        logw("Counting finished. Starting fit...", Some("info"));
        let mut fit = SpectrumFitter::new();
        let result = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec());
        if result.is_ok() {
            minc = result.unwrap() as u16;
        } else {
            logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
            minc = 5;
        }

        // minc = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;

        if minc <= 0 {
            panic!("Fitted min_count value is zero or negative!");
        }

        logw(format!("Fit done! Fitted min_count value: {}. Starting filtering...", minc).as_str(), Some("info"));

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            false
        });
    } else {
        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }


            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize- 1].saturating_add(1);
            }

            false
        });
    }

    outdict.shrink_to_fit();
    minmaxdict.shrink_to_fit();
    histovec.shrink_to_fit();
    drop(countmap);

    (outdict, minmaxdict, themap, histovec, minc)
}


#[cfg(feature = "wasm")]
fn bloom_filter_preprocessing_wasm<IntT>(
    file1:    &mut WebSysFile,
    file2:    &mut WebSysFile,
    k:         usize,
    qual:     &QualOpts,
    do_fit:    bool,
) -> (HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u32>, u16)
where
    IntT: for<'a> UInt<'a>,
{
    logw("Getting kmers from first file with Bloom filter. Creating reader and initialising filter...", Some("info"));

    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    let mut themap     = HashMap::with_hasher(BuildHasherDefault::default());

    let mut reader = open_fastq(file1);
    let mut histovec : Vec<u32> = vec![0; MAXSIZEHISTO];

    let mut kmer_filter = KmerFilter::new(qual.min_count);
    kmer_filter.init();

    logw("Entering while loop for the first file...", Some("info"));

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer(); // TODO: potential really small improvement, get only hash for the bloom filter, then get the rest.
            if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                    outdict.entry(hc).or_insert(km);
                    minmaxdict.entry(hnc).or_insert(hc);
                }
            }
        }


    }
    logw("Finished getting kmers from first file. Starting with the second...", Some("info"));

    let mut reader = open_fastq(file2);

    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        // put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                    outdict.entry(hc).or_insert(km);
                    minmaxdict.entry(hnc).or_insert(hc);
                }
            }
        }
    }

    logw("Finished getting kmers from the second file", Some("info"));
    logw("Second part of filtering...", Some("info"));

    // // Now, get themap, histovec, and filter outdict and minmaxdict
    let countmap = kmer_filter.get_counts_map();
    countmap.shrink_to_fit();


    // countmap.retain(|h, tup| {
    //     if tup.0 >= qual.min_count {
    //         themap
    //             .entry(*h)
    //             .or_insert( RefCell::new(HashInfoSimple {
    //                         hnc:    tup.1,
    //                         b:      tup.2,
    //                         pre:    Vec::new(),
    //                         post:   Vec::new(),
    //                         counts: tup.0,
    //                 }) );
    //     } else {
    //         outdict.remove(h);
    //         minmaxdict.remove(&tup.1);
    //     }
    //
    //     if tup.0 as usize > MAXSIZEHISTO {
    //         histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
    //     } else {
    //         histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
    //     }
    //
    //     false
    // });
    //
    // outdict.shrink_to_fit();
    // minmaxdict.shrink_to_fit();

    // coses

    // This can be optimised. also better written: I had to repeat the code for the retains, to try to improve slightly the running time in
    // case no autofitting is requested. In any case, it could be improved in the future.
    let mut minc = qual.min_count;
    if do_fit {
        for (_, tup) in countmap.iter() {
            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
        }

        // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
        // large (and so that we can detect it).
        logw("Counting finished. Starting fit...", Some("info"));
        let mut fit = SpectrumFitter::new();
        let result = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec());
        if result.is_ok() {
            minc = result.unwrap() as u16;
        } else {
            logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
            minc = 5;
        }

        // minc = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;

        if minc <= 0 {
            panic!("Fitted min_count value is zero or negative!");
        }

        logw(format!("Fit done! Fitted min_count value: {}. Starting filtering...", minc).as_str(), Some("info"));

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            false
        });
    } else {
        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
            false
        });
    }

    outdict.shrink_to_fit();
    minmaxdict.shrink_to_fit();
    // coses


    (outdict, minmaxdict, themap, histovec, minc)
}


#[cfg(not(feature = "wasm"))]
fn bloom_filter_preprocessing_standalone<IntT>(
    file1:  &str,
    file2:  &str,
    k:      usize,
    qual:   &QualOpts,
    do_fit: bool,
    out_path: &mut Option<PathBuf>,
) -> (HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    log::info!("Getting kmers from first file with Bloom filter. Creating reader and initialising filter...");

    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    let mut themap     = HashMap::with_hasher(BuildHasherDefault::default());

    let mut reader =
        parse_fastx_file(file1).unwrap_or_else(|_| panic!("Invalid path/file: {file1}"));
    let mut histovec : Vec<u32> = vec![0; MAXSIZEHISTO];

    let mut kmer_filter = KmerFilter::new(qual.min_count);
    kmer_filter.init();

    logw(format!("Entering while loop for the first file...").as_str(), Some("info"));

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer(); // TODO: potential really small improvement, get only hash for the bloom filter, then get the rest.
            if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                    outdict.entry(hc).or_insert(km);
                    minmaxdict.entry(hnc).or_insert(hc);
                }
            }
        }


    }
    log::info!("Finished getting kmers from first file. Starting with the second...");

    reader =
        parse_fastx_file(file2).unwrap_or_else(|_| panic!("Invalid path/file: {file2}"));

    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        // put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                if Ordering::is_eq(kmer_filter.filter(hc, hnc, b)) {
                    outdict.entry(hc).or_insert(km);
                    minmaxdict.entry(hnc).or_insert(hc);
                }
            }
        }
    }

    log::info!("Finished getting kmers from the second file");
    log::info!("Second part of filtering...");

    // Now, get themap, histovec, and filter outdict and minmaxdict
    let countmap = kmer_filter.get_counts_map();
    countmap.shrink_to_fit();
    let minc;

    // This can be optimised. also better written: I had to repeat the code for the retains, to try to improve slightly the running time in
    // case no autofitting is requested. In any case, it could be improved in the future.
    if do_fit {
        for (_, tup) in countmap.iter() {
            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
        }

        // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
        // large (and so that we can detect it).
        log::info!("Counting finished. Starting fit...");
        let mut fit = SpectrumFitter::new();

        let result = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec());
        if result.is_ok() {
            minc = result.unwrap() as u16;
        } else {
            logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
            minc = 5;
        }

        // let minc = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;

        if minc <= 0 {
            panic!("Fitted min_count value is zero or negative!");
        }

        log::info!("Fit done! Fitted min_count value: {}. Starting filtering...", minc);

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            false
        });
    } else {
        minc = qual.min_count;

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
            false
        });
    }

    outdict.shrink_to_fit();
    minmaxdict.shrink_to_fit();

    if out_path.is_some() {
        // Plotting!
        let backend = BitMapBackend::new(out_path.as_ref().unwrap().as_path(), (1280, 960));

        let root = backend.into_drawing_area();

        let _ = root.fill(&WHITE);

        let mut chart = ChartBuilder::on(&root)
            .x_label_area_size(35)
            .y_label_area_size(40)
            .margin(5)
            // .caption("k-mer spectrum", ("sans-serif", 30.0))
            .caption("k-mer spectrum", ("ibm-plex-sans", 30.0))
            .build_cartesian_2d((0u32..(MAXSIZEHISTO as u32)).into_segmented(), 0u32..200000u32).unwrap();

        chart
            .configure_mesh()
            .disable_x_mesh()
            .bold_line_style(WHITE.mix(0.3))
            .y_desc("Counts")
            .x_desc("k-mer frequency")
            // .axis_desc_style(("sans-serif", 15))
            .axis_desc_style(("ibm-plex-sans", 15))
            .draw().unwrap();

        chart.draw_series(
            Histogram::vertical(&chart)
                .style(RED.filled())
                .data(histovec.iter().enumerate().map(|(i, x)| (i as u32, *x))),
        ).unwrap();



        // TODO: I have spent too much time trying to draw a vertical line or a rectangle to show in the
        //       histogram the fitted limit. Not more at least until I decide to lose more of my life.
        // https://stackoverflow.com/questions/78776201/how-to-dynamically-use-plotter-segmentvalue

        // backend.draw_rect(
        //     (0i32, 200000i32),
        //     (fitted_min_count as i32, 0i32),
        //     &BLACK,
        //     true,
        // ).unwrap();

        // let testnum = fitted_min_count as i32;
        // let rectangle = Rectangle::new(
        //     [(0, 200000), (5, 0)],
        //     BLUE.mix(0.5).filled(),
        // );
        //
        // chart.draw_series(std::iter::once(rectangle.into_dyn())).unwrap();

        // chart.plotting_area().draw(&rectangle).unwrap();

        // chart.draw_series(LineSeries::new(
        //     [(0i32, 0i32), (fitted_min_count as i32, 200000i32)].iter(),
        //     &BLUE,
        // )).unwrap();

        root.present().expect("Unable to write result to file. Does the output folder exist?");
    }

    (outdict, minmaxdict, themap)
}


#[cfg(not(feature = "wasm"))]
fn get_map_with_counts(
    invec:      &Vec<(u64, u64, u8)>,
    min_count:  u16,
    out_path:   &mut Option<PathBuf>,
) -> HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>
{
    let mut outdict = HashMap::with_hasher(BuildHasherDefault::default());

    let mut i = 0;
    let mut c : u16 = 0;
    let mut tmphash = invec[i].0;
    // let mut tmpcounter = 0;

    // I'm not sure if there is another way of avoiding the extra conditional checks for plotting (or perhaps the
    // compiler is smart enough to write in assembly exactly what I am now writing here?) than to rewrite the function
    // with a big "if" as it is here, to gain a bit of efficiency (not sure how much, but well).
    if out_path.is_some() {
        let mut plotvec : Vec<u16> = Vec::new(); // For plotting

        while i < invec.len() {
            if tmphash != invec[i].0 {
                if c >= min_count {
                    // tmpcounter += 1;
                    outdict
                        .entry(tmphash)
                        .or_insert( RefCell::new(HashInfoSimple {
                                    hnc:    invec[i - 1].1,
                                    b:      invec[i - 1].2,
                                    pre:    Vec::new(),
                                    post:   Vec::new(),
                                    counts: c,
                            }) );
                } else {
                    plotvec.push(c);
                }
                tmphash = invec[i].0;
                c = 1;
            } else {
                c = c.saturating_add(1);
            }
            i += 1;
        }


        if c >= min_count {
            // tmpcounter += 1;
            outdict
                .entry(tmphash)
                .or_insert( RefCell::new(HashInfoSimple {
                            hnc:    invec[i - 1].1,
                            b:      invec[i - 1].2,
                            pre:    Vec::new(),
                            post:   Vec::new(),
                            counts: c,
                    }) );
        } else {
            plotvec.push(c);
        }
        plotvec.shrink_to_fit();

        // Plotting!

        let root = BitMapBackend::new(out_path.as_ref().unwrap().as_path(), (1280, 960)).into_drawing_area();

        let _ = root.fill(&WHITE);

        let mut chart = ChartBuilder::on(&root)
            .x_label_area_size(35)
            .y_label_area_size(40)
            .margin(5)
            .caption("k-mer spectrum", ("sans-serif", 30.0))
            .build_cartesian_2d((0u32..200u32).into_segmented(), 0u32..200000u32).unwrap();

        chart
            .configure_mesh()
            .disable_x_mesh()
            .bold_line_style(WHITE.mix(0.3))
            .y_desc("Counts")
            .x_desc("k-mer frequency")
            .axis_desc_style(("sans-serif", 15))
            .draw().unwrap();

        chart.draw_series(
            Histogram::vertical(&chart)
                .style(RED.filled())
                .data(plotvec.iter().map(|x: &u16| (*x as u32, 1)).chain(outdict.iter().map(|(_, x)| (x.borrow().counts as u32, 1)))),
        ).unwrap();

        root.present().expect("Unable to write result to file. Does the output folder exist?");

    //     exit(1);

        // log::debug!("Good kmers {}", tmpcounter);
    } else {
        // Here we don't need to check for plotting or anything
        while i < invec.len() {
            if tmphash != invec[i].0 {
                if c >= min_count {
                    // tmpcounter += 1;
                    outdict
                        .entry(tmphash)
                        .or_insert( RefCell::new(HashInfoSimple {
                                    hnc:    invec[i - 1].1,
                                    b:      invec[i - 1].2,
                                    pre:    Vec::new(),
                                    post:   Vec::new(),
                                    counts: c,
                            }) );
                }
                tmphash = invec[i].0;
                c = 1;
            } else {
                c = c.saturating_add(1);
            }
            i += 1;
        }


        if c >= min_count {
            // tmpcounter += 1;
            outdict
                .entry(tmphash)
                .or_insert( RefCell::new(HashInfoSimple {
                            hnc:    invec[i - 1].1,
                            b:      invec[i - 1].2,
                            pre:    Vec::new(),
                            post:   Vec::new(),
                            counts: c,
                    }) );
        }
    }
    return outdict;
}



#[cfg(not(feature = "wasm"))]
fn get_map_with_counts_and_fit(
    invec:      &mut Vec<(u64, u64, u8)>,
    out_path:   &mut Option<PathBuf>,
) -> HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>
{
    let mut outdict = HashMap::with_hasher(BuildHasherDefault::default());

    let mut i = 0;
    let mut c : u16 = 0;
    let mut tmphash = invec[i].0;
    // let mut tmpcounter = 0;

    // Here we don't have the same issue as with the pre-defined min_count setting.
    let mut plotvec : Vec<u32> = vec![0 as u32; MAXSIZEHISTO]; // For plotting

    // We need to construct the histogram as well

    while i < invec.len() {
        if tmphash != invec[i].0 {
            // tmpcounter += 1;
            outdict
                .entry(tmphash)
                .or_insert( RefCell::new(HashInfoSimple {
                            hnc:    invec[i - 1].1,
                            b:      invec[i - 1].2,
                            pre:    Vec::new(),
                            post:   Vec::new(),
                            counts: c,
                    }) );

            if c as usize > MAXSIZEHISTO {
                plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
            }

            tmphash = invec[i].0;
            c = 1;

        } else {
            c = c.saturating_add(1);
        }
        i += 1;
    }


    // tmpcounter += 1;
    outdict
        .entry(tmphash)
        .or_insert( RefCell::new(HashInfoSimple {
                    hnc:    invec[i - 1].1,
                    b:      invec[i - 1].2,
                    pre:    Vec::new(),
                    post:   Vec::new(),
                    counts: c,
            }) );


    if c as usize > MAXSIZEHISTO {
        plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
    } else {
        plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
    }

    invec.clear(); invec.shrink_to_fit(); // Quick optimisation

    // We need to do the fit!
    log::info!("Counting finished. Starting fit...");
    let mut fit = SpectrumFitter::new();
    // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
    // large (and so that we can detect it).
    // let fitted_min_count = fit.fit_histogram(plotvec[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;
    let fitted_min_count : u16;

    let result = fit.fit_histogram(plotvec[..(MAXSIZEHISTO - 1)].to_vec());
    if result.is_ok() {
        fitted_min_count = result.unwrap() as u16;
    } else {
        logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
        fitted_min_count = 5;
    }

    if fitted_min_count <= 0 {
        panic!("Fitted min_count value is zero or negative!");
    }

    log::info!("Fit done! Fitted min_count value: {}. Starting filtering...", fitted_min_count);

    outdict.retain(|_, rc| {
        rc.borrow().counts >= fitted_min_count
    });
    outdict.shrink_to_fit();

    if out_path.is_some() {
        // Plotting!
        let backend = BitMapBackend::new(out_path.as_ref().unwrap().as_path(), (1280, 960));

        let root = backend.into_drawing_area();

        let _ = root.fill(&WHITE);

        let mut chart = ChartBuilder::on(&root)
            .x_label_area_size(35)
            .y_label_area_size(40)
            .margin(5)
            // .caption("k-mer spectrum", ("sans-serif", 30.0))
            .caption("k-mer spectrum", ("ibm-plex-sans", 30.0))
            .build_cartesian_2d((0u32..(MAXSIZEHISTO as u32)).into_segmented(), 0u32..200000u32).unwrap();

        chart
            .configure_mesh()
            .disable_x_mesh()
            .bold_line_style(WHITE.mix(0.3))
            .y_desc("Counts")
            .x_desc("k-mer frequency")
            // .axis_desc_style(("sans-serif", 15))
            .axis_desc_style(("ibm-plex-sans", 15))
            .draw().unwrap();

        chart.draw_series(
            Histogram::vertical(&chart)
                .style(RED.filled())
                .data(plotvec.iter().enumerate().map(|(i, x)| (i as u32, *x))),
        ).unwrap();



        // TODO: I have spent too much time trying to draw a vertical line or a rectangle to show in the
        //       histogram the fitted limit. Not more at least until I decide to lose more of my life.
        // https://stackoverflow.com/questions/78776201/how-to-dynamically-use-plotter-segmentvalue

        // backend.draw_rect(
        //     (0i32, 200000i32),
        //     (fitted_min_count as i32, 0i32),
        //     &BLACK,
        //     true,
        // ).unwrap();

        // let testnum = fitted_min_count as i32;
        // let rectangle = Rectangle::new(
        //     [(0, 200000), (5, 0)],
        //     BLUE.mix(0.5).filled(),
        // );
        //
        // chart.draw_series(std::iter::once(rectangle.into_dyn())).unwrap();

        // chart.plotting_area().draw(&rectangle).unwrap();

        // chart.draw_series(LineSeries::new(
        //     [(0i32, 0i32), (fitted_min_count as i32, 200000i32)].iter(),
        //     &BLUE,
        // )).unwrap();

        root.present().expect("Unable to write result to file. Does the output folder exist?");
    }


    // log::debug!("Good kmers {}", tmpcounter);
//     exit(1);

    outdict
}



#[cfg(feature = "wasm")]
fn get_map_wasm(
    invec:      &mut Vec<(u64, u64, u8)>,
    min_count:  u16,
    do_fit:     bool,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u32>, u16)
{
    let mut outdict = HashMap::with_hasher(BuildHasherDefault::default());


    let mut i = 0;
    let mut c : u16 = 0;
    let mut minc : u16 = min_count;
    let mut tmphash = invec[i].0;
    // let mut tmpcounter = 0;

    let mut plotvec : Vec<u32> = vec![0 as u32; MAXSIZEHISTO]; // For plotting

    if do_fit {
        while i < invec.len() {
            if tmphash != invec[i].0 {
                // tmpcounter += 1;
                outdict
                    .entry(tmphash)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    invec[i - 1].1,
                                b:      invec[i - 1].2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: c,
                        }) );

                if c as usize > MAXSIZEHISTO {
                    plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
                } else {
                    plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
                }

                tmphash = invec[i].0;
                c = 1;

            } else {
                c = c.saturating_add(1);
            }
            i += 1;
        }


        // tmpcounter += 1;
        outdict
            .entry(tmphash)
            .or_insert( RefCell::new(HashInfoSimple {
                        hnc:    invec[i - 1].1,
                        b:      invec[i - 1].2,
                        pre:    Vec::new(),
                        post:   Vec::new(),
                        counts: c,
                }) );


        if c as usize > MAXSIZEHISTO {
            plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
        } else {
            plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
        }

        invec.clear(); invec.shrink_to_fit(); // Quick optimisation

        // We need to do the fit!
        logw("Counting finished. Starting fit...", Some("info"));
        let mut fit = SpectrumFitter::new();
        // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
        // large (and so that we can detect it).
        // minc = fit.fit_histogram(plotvec[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;

        let result = fit.fit_histogram(plotvec[..(MAXSIZEHISTO - 1)].to_vec());
        if result.is_ok() {
            minc = result.unwrap() as u16;
        } else {
            logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
            minc = 5;
        }

        if minc <= 0 {
            panic!("Fitted min_count value is zero or negative!");
        }

        logw(format!("Fit done! Fitted min_count value: {}. Starting filtering...", minc).as_str(), Some("info"));

        outdict.retain(|_, rc| {
            rc.borrow().counts >= minc
        });
        outdict.shrink_to_fit();

    } else {
        while i < invec.len() {
            if tmphash != invec[i].0 {
                if c >= minc {
                    // tmpcounter += 1;
                    outdict
                        .entry(tmphash)
                        .or_insert( RefCell::new(HashInfoSimple {
                                    hnc:    invec[i - 1].1,
                                    b:      invec[i - 1].2,
                                    pre:    Vec::new(),
                                    post:   Vec::new(),
                                    counts: c,
                            }) );
                }

                if c as usize > MAXSIZEHISTO {
                    plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
                } else {
                    plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
                }

                tmphash = invec[i].0;
                c = 1;
            } else {
                c = c.saturating_add(1);
            }
            i += 1;
        }


        if c >= minc {
            // tmpcounter += 1;
            outdict
                .entry(tmphash)
                .or_insert( RefCell::new(HashInfoSimple {
                            hnc:    invec[i - 1].1,
                            b:      invec[i - 1].2,
                            pre:    Vec::new(),
                            post:   Vec::new(),
                            counts: c,
                    }) );
        }

        if c as usize > MAXSIZEHISTO {
            plotvec[MAXSIZEHISTO - 1] = plotvec[MAXSIZEHISTO - 1].saturating_add(1);
        } else {
            plotvec[c as usize - 1] = plotvec[c as usize - 1].saturating_add(1);
        }
        // logw(format!("Good kmers {}", tmpcounter).as_str(), Some("debug"));
    }
    return (outdict, plotvec, minc);
}


fn update_countmap(
    invec    : &    Vec<(u64, u64, u8)>,
    countmap : &mut HashMap::<u64, (u16, u64, u8), BuildHasherDefault<NoHashHasher<u64>>>,
) {

    let mut i = 0;
    let mut c : u16 = 0;
    let mut tmphash = invec[i].0;
    // let mut tmpcounter = 0;


    while i < invec.len() {
        if tmphash != invec[i].0 {

            let tmpref = countmap.entry(tmphash).or_insert( (0, invec[i - 1].1, invec[i - 1].2) );
            tmpref.0 = tmpref.0.saturating_add(c);

            tmphash = invec[i].0;
            c = 1;
        } else {
            c = c.saturating_add(1);
        }
        i += 1;
    }

    let tmpref = countmap.entry(tmphash).or_insert( (0, invec[i - 1].1, invec[i - 1].2) );
    tmpref.0 = tmpref.0.saturating_add(c);
}


#[cfg(not(feature = "wasm"))]
fn chunked_processing_standalone<IntT>(
    file1:    &str,
    file2:    &str,
    k:        usize,
    qual:     &QualOpts,
    outvec:   &mut Vec<(u64, u64, u8)>,
    csize:    usize,
    do_fit:   bool,
    out_path: &mut Option<PathBuf>,
) -> (HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    log::info!("Getting kmers from first file. Creating reader...");

    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    let mut themap     = HashMap::with_hasher(BuildHasherDefault::default());
    let mut countmap : HashMap::<u64, (u16, u64, u8), BuildHasherDefault<NoHashHasher<u64>>> = HashMap::with_hasher(BuildHasherDefault::default());

    log::info!("Getting kmers from file {file1}. Creating reader...");
    let mut reader =
        parse_fastx_file(file1).unwrap_or_else(|_| panic!("Invalid path/file: {file1}"));

    let mut histovec : Vec<u32> = vec![0; MAXSIZEHISTO];

    log::info!("Entering while loop...");

    let mut i_record = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
        }

        i_record += 1;
        if i_record >= csize {
            // Processssssss! And reset.
            if !outvec.is_empty() {
                log::info!("Processing chunk. Sorting k-mers...");
                outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                log::info!("k-mers sorted. Counting k-mers...");
                // Then, do a counting of everything and save the results in a dictionary and return it

                update_countmap(&outvec, &mut countmap);
            }

            // Reset
            outvec.clear();
            i_record = 0;
        }

    }
    log::info!("Finished getting kmers from first file. Starting with the second...");

    reader =
        parse_fastx_file(file2).unwrap_or_else(|_| panic!("Invalid path/file: {file2}"));

    // Filling the seq of the second file!
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        // put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            seqrec.seq(),
            rl,
            seqrec.qual(),
            k,
            qual.min_qual,
            true,
        );
        if let Some(mut kmer_it) = kmer_opt {
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            outdict.entry(hc).or_insert(km);
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                outdict.entry(hc).or_insert(km);
                minmaxdict.entry(hnc).or_insert(hc);
            }
        }

        i_record += 1;
        if i_record >= csize {
            // Processssssss! And reset.
            if !outvec.is_empty() {
                log::info!("Processing chunk. Sorting k-mers...");
                outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                log::info!("k-mers sorted. Counting k-mers...");
                // Then, do a counting of everything and save the results in a dictionary and return it

                update_countmap(&outvec, &mut countmap);
            }

            // Reset
            outvec.clear();
            i_record = 0;
        }
    }

    if i_record > 0 {
        // Processssssss! And reset.
        if !outvec.is_empty() {
            log::info!("Processing last chunk. Sorting k-mers...");
            outvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            log::info!("k-mers sorted. Counting k-mers...");
            // Then, do a counting of everything and save the results in a dictionary and return it

            update_countmap(&outvec, &mut countmap);
        }
        // Reset
        outvec.clear();
    }
    log::info!("Finished getting kmers from the second file");
    log::info!("Filtering...");

    // Now, get themap, histovec, and filter outdict and minmaxdict
    countmap.shrink_to_fit();
    let minc;

    // This can be optimised. also better written: I had to repeat the code for the retains, to try to improve slightly the running time in
    // case no autofitting is requested. In any case, it could be improved in the future.
    if do_fit {
        for (_, tup) in countmap.iter() {
            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
        }
        // Remove the last bin, as it might affect the fit, but we want it in the vector to plot it in case the coverage is really
        // large (and so that we can detect it).
        log::info!("Counting finished. Starting fit...");
        let mut fit = SpectrumFitter::new();
        // let minc = fit.fit_histogram(histovec.clone()[..(MAXSIZEHISTO - 1)].to_vec()).expect("Fit to the k-mer spectrum failed!") as u16;

        let result = fit.fit_histogram(histovec[..(MAXSIZEHISTO - 1)].to_vec());
        if result.is_ok() {
            minc = result.unwrap() as u16;
        } else {
            logw("Fit has not converged. The default value of 5 will be used. You should check whether this value is appropiated or not by looking at the k-mer spectrum histogram.", Some("warn"));
            minc = 5;
        }

        if minc <= 0 {
            panic!("Fitted min_count value is zero or negative!");
        }

        log::info!("Fit done! Fitted min_count value: {}. Starting filtering...", minc);

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            false
        });
    } else {
        minc = qual.min_count;

        countmap.retain(|h, tup| {
            if tup.0 >= minc {
                themap
                    .entry(*h)
                    .or_insert( RefCell::new(HashInfoSimple {
                                hnc:    tup.1,
                                b:      tup.2,
                                pre:    Vec::new(),
                                post:   Vec::new(),
                                counts: tup.0,
                        }) );
            } else {
                outdict.remove(h);
                minmaxdict.remove(&tup.1);
            }

            if tup.0 as usize > MAXSIZEHISTO {
                histovec[MAXSIZEHISTO - 1] = histovec[MAXSIZEHISTO - 1].saturating_add(1);
            } else {
                histovec[tup.0 as usize - 1] = histovec[tup.0 as usize - 1].saturating_add(1);
            }
            false
        });
    }

    drop(countmap);
    outdict.shrink_to_fit();
    minmaxdict.shrink_to_fit();

    if out_path.is_some() {
        // Plotting!
        let backend = BitMapBackend::new(out_path.as_ref().unwrap().as_path(), (1280, 960));

        let root = backend.into_drawing_area();

        let _ = root.fill(&WHITE);

        let mut chart = ChartBuilder::on(&root)
            .x_label_area_size(35)
            .y_label_area_size(40)
            .margin(5)
            // .caption("k-mer spectrum", ("sans-serif", 30.0))
            .caption("k-mer spectrum", ("ibm-plex-sans", 30.0))
            .build_cartesian_2d((0u32..(MAXSIZEHISTO as u32)).into_segmented(), 0u32..200000u32).unwrap();

        chart
            .configure_mesh()
            .disable_x_mesh()
            .bold_line_style(WHITE.mix(0.3))
            .y_desc("Counts")
            .x_desc("k-mer frequency")
            // .axis_desc_style(("sans-serif", 15))
            .axis_desc_style(("ibm-plex-sans", 15))
            .draw().unwrap();

        chart.draw_series(
            Histogram::vertical(&chart)
                .style(RED.filled())
                .data(histovec.iter().enumerate().map(|(i, x)| (i as u32, *x))),
        ).unwrap();



        // TODO: I have spent too much time trying to draw a vertical line or a rectangle to show in the
        //       histogram the fitted limit. Not more at least until I decide to lose more of my life.
        // https://stackoverflow.com/questions/78776201/how-to-dynamically-use-plotter-segmentvalue

        // backend.draw_rect(
        //     (0i32, 200000i32),
        //     (fitted_min_count as i32, 0i32),
        //     &BLACK,
        //     true,
        // ).unwrap();

        // let testnum = fitted_min_count as i32;
        // let rectangle = Rectangle::new(
        //     [(0, 200000), (5, 0)],
        //     BLUE.mix(0.5).filled(),
        // );
        //
        // chart.draw_series(std::iter::once(rectangle.into_dyn())).unwrap();

        // chart.plotting_area().draw(&rectangle).unwrap();

        // chart.draw_series(LineSeries::new(
        //     [(0i32, 0i32), (fitted_min_count as i32, 200000i32)].iter(),
        //     &BLUE,
        // )).unwrap();

        root.present().expect("Unable to write result to file. Does the output folder exist?");
    }

    (outdict, minmaxdict, themap)
}


/// Read fastq files, get the reads, get the k-mers, count them, filter them by count, and get some way of recovering the sequence later.
#[cfg(not(feature = "wasm"))]
pub fn preprocessing_standalone<IntT>(
    input_files:    &[InputFastx],
    k:              usize,
    qual:           &QualOpts,
    timevec:        &mut Vec<Instant>,
    out_path:       &mut Option<PathBuf>,
    csize   :       usize,
    do_bloom:       bool,
    do_fit  :       bool,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u64>, HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    log::info!("Starting preprocessing_standalone with k = {k}");

    if do_bloom {
        // Build indexes
        log::info!("Processing using a Bloom filter");

        let (thedict, maxmindict, themap) = bloom_filter_preprocessing_standalone::<IntT>(
            &input_files[0].1,
            input_files[0].2.as_ref().expect("No paired reads!"),
            k,
            qual,
            do_fit,
            out_path,
        );
        return (themap, Vec::new(), thedict, maxmindict);

    } else if csize <= 0 {
        if any_fastq(input_files) {
            log::info!("FASTQ files filtered with: {qual}");
        } else {
            panic!("Input files are not FASTQ");
        }

        let total_size  = input_files.len();
        if total_size > 1 {panic!("Not expecting more than a pair of pair-end reads right now!");};


        // First, we want to fill our mega-vector with all k-mers from both paired-end reads
        log::info!("Filling vector");
        println!("checks enabled!");

        let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
        let (theseq, thedict, maxmindict) = get_kmers_from_both_files_and_the_dict_and_the_seq::<IntT>(&input_files[0].1,
            input_files[0].2.as_ref().expect("No paired reads!"),
            k,
            qual,
            &mut tmpvec);

        exit(0);
        timevec.push(Instant::now());
        log::info!("Kmers extracted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());


        // Then, we want to sort it according to the hash
        log::debug!("Number of kmers BEFORE cleaning: {:?}", tmpvec.len());
        log::info!("Sorting vector");
        tmpvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        timevec.push(Instant::now());
        log::info!("Kmers sorted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());

        // Then, do a counting of everything and save the results in a dictionary and return it
        timevec.push(Instant::now());
        log::info!("Counting k-mers");
        let themap;

        if !do_fit {
            themap = get_map_with_counts(&mut tmpvec, qual.min_count, out_path);
        } else {
            themap = get_map_with_counts_and_fit(&mut tmpvec, out_path);
        }
        drop(tmpvec);

        timevec.push(Instant::now());
        log::info!("Kmers counted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());

        return (themap, theseq, thedict, maxmindict)
    } else {
        // Build indexes
        log::info!("Processing in chunks of size {} the input files", csize);

        let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
        let (thedict, maxmindict, themap) = chunked_processing_standalone::<IntT>(
            &input_files[0].1,
            input_files[0].2.as_ref().expect("No paired reads!"),
            k,
            qual,
            &mut tmpvec,
            csize,
            do_fit,
            out_path,
        );
        drop(tmpvec);
        return (themap, Vec::new(), thedict, maxmindict);
    }

}



#[cfg(feature = "wasm")]
/// Main preprocessing function for wasm
pub fn preprocessing_wasm<IntT>(
    file1   : &mut WebSysFile,
    file2   : &mut WebSysFile,
    k       : usize,
    qual    : &QualOpts,
    csize   : usize,
    do_bloom: bool,
    do_fit  : bool,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Option<HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u32>, u16)
where
    IntT: for<'a> UInt<'a>,
{
    if do_bloom {
        // Build indexes
        logw("Processing using a Bloom filter", Some("info"));

        let (thedict, maxmindict, themap, histovec, used_min_count) = bloom_filter_preprocessing_wasm::<IntT>(
            file1,
            file2,
            k,
            qual,
            do_fit
        );
        return (themap, Some(thedict), maxmindict, histovec, used_min_count);

    } else if csize <= 0 {
        // Build indexes
        logw("Starting preprocessing with k = {k}", Some("info"));

        // First, we want to fill our mega-vector with all k-mers from both paired-end reads
        logw("Filling vector", Some("info"));

        let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
        let (thedict, maxmindict) = get_kmers_from_both_files_wasm::<IntT>(file1,
                                                                           file2,
                                                                           k,
                                                                           qual,
                                                                           &mut tmpvec
        );

        logw("Kmers extracted", Some("info"));


        // Then, we want to sort it according to the hash
        // log::debug!("Number of kmers BEFORE cleaning: {:?}", tmpvec.len());
        logw("Sorting vector", Some("info"));
        tmpvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        logw("Kmers sorted.", Some("info"));


        // Then, do a counting of everything and save the results in a dictionary and return it
        logw("Counting k-mers", Some("info"));
        let (themap, mut histovec, used_min_count) = get_map_wasm(&mut tmpvec, qual.min_count, do_fit);
        histovec.shrink_to_fit();
        drop(tmpvec);

        logw("Kmers counted.", Some("info"));

        return (themap, Some(thedict), maxmindict, histovec, used_min_count);
    } else {
        // Build indexes
        logw(format!("Processing in chunks of size {} the input files", csize).as_str(), Some("info"));

        let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
        let (thedict, maxmindict, themap, mut histovec, used_min_count) = chunked_processing_wasm::<IntT>(
            file1,
            file2,
            k,
            qual,
            &mut tmpvec,
            csize,
            do_fit,
        );
        drop(tmpvec);
        histovec.shrink_to_fit();
        return (themap, Some(thedict), maxmindict, histovec, used_min_count);
    }
}
