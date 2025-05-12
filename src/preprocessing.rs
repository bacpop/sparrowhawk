//! Some docs should be here

use core::panic;
use std::time::Instant;

use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault, cell::*};

use rayon::prelude::*;

#[cfg(not(feature = "wasm"))]
use indicatif::ProgressIterator;

#[cfg(not(feature = "wasm"))]
use needletail::parse_fastx_file;

// use std::process::exit;

#[cfg(not(feature = "wasm"))]
use plotters::prelude::*;

use project_root;

use super::QualOpts;
use super::HashInfoSimple;
use super::bit_encoding::{encode_base, rc_base};
use crate::io_utils::any_fastq;
use crate::bit_encoding::UInt;
use crate::kmer::Kmer;
use crate::loG;


/// Tuple for name and fasta or paired fastq input
pub type InputFastx = (String, String, Option<String>);

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;
#[cfg(feature = "wasm")]
use wasm_bindgen_file_reader::WebSysFile;
#[cfg(feature = "wasm")]
use crate::fastx_wasm::open_fastq;
#[cfg(feature = "wasm")]
use seq_io::fastq::Record;

// =====================================================================================================


fn put_these_nts_into_an_efficient_vector(charseq : &[u8], compseq : &mut Vec<u64>, occ : u8) {
    let mut tmpu64 : u64 = 0;
    let mut tmpind : u8  = 0;

    if occ != 0 {
        tmpind = occ;
        tmpu64 = compseq.pop().unwrap();
    }
//     println!("{}", tmpind);

    for nt in charseq {
//         println!("\n{:#010b}\n{}\n{:#066b}\n{:#066b}", *nt, tmpind, tmpu64, (encode_base(*nt) as u64));
        tmpu64 <<= 2;
        tmpu64 |= encode_base(*nt) as u64;
//         println!("{:#066b}", tmpu64);
        if tmpind == 31 {
            compseq.push(tmpu64);
            tmpu64 = 0;
            tmpind = 0;
        } else {
            tmpind += 1;
        }
//         println!("{}", tmpind);
    }

    if tmpind != 0 {
        compseq.push(tmpu64);
    }
}

fn put_these_nts_into_an_efficient_vector_rc(charseq : &[u8], compseq : &mut Vec<u64>, occ : u8) {
    let mut tmpu64 : u64 = 0;
    let mut tmpind : u8  = 0;

    if occ != 0 {
        tmpind = occ;
        tmpu64 = compseq.pop().unwrap();
    }
//     println!("{}", tmpind);
    for nt in charseq.iter().rev() {
//         println!("\n{:#010b}\n{}\n{:#066b}\n{:#066b}", *nt, tmpind, tmpu64, (rc_base(encode_base(*nt)) as u64));
        tmpu64 <<= 2;
        tmpu64 |= rc_base(encode_base(*nt)) as u64;
//         println!("{:#066b}", tmpu64);
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
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                println!("\n\t- COLLISIONS 1 !!! Hash: {:?}", hc);
                println!("{:#0258b}\n{:#0258b}", *testkm, km);
            }
            // } else {
            //     println!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    println!("\n\t- COLLISIONS 2 !!! Hash: {:?}", hc);
                    println!("{:#0258b}\n{:#0258b}", *testkm, km);
                }
                // } else {
                //     println!("\n\t\t- NOT COLLISIONS!!!");
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

//     println!("MITAD Length of seq. vec.: {}, total length of first files: {}, THING {}, number of recs: {}", theseq.len(), itrecord, (itrecord % 32) as u8, realit);
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
//         println!("{}", seqrec.num_bases());
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
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                println!("\n\t- COLLISIONS 3 !!! Hash: {:?}", hc);
                println!("{:#0258b}\n{:#0258b}", *testkm, km);
            }
            // } else {
            //     println!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    println!("\n\t- COLLISIONS 4 !!! Hash: {:?}", hc);
                    println!("{:#0258b}\n{:#0258b}", *testkm, km);
                }
                // } else {
                //     println!("\n\t\t- NOT COLLISIONS!!!");
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
//         println!("{:#066b}", tmpu64);
        theseq.push(tmpu64);
    }

    log::info!("Finished getting kmers from file {filename2}");
    println!("Length of seq. vec.: {}, total length of both files: {}", theseq.len(), itrecord);


    (theseq, outdict, minmaxdict)
}



#[cfg(feature = "wasm")]
fn get_kmers_from_both_files_wasm<IntT>(
    file1:    &mut WebSysFile,
    file2:    &mut WebSysFile,
    k:         usize,
    qual:      &QualOpts,
    outvec:    &mut Vec<(u64, u64, u8)>,
) -> (Vec<u64>, HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    let mut outdict    = HashMap::with_hasher(BuildHasherDefault::default());
    let mut minmaxdict = HashMap::with_hasher(BuildHasherDefault::default());
    loG("Getting kmers from first file. Creating reader...", Some("info"));
    let mut reader = open_fastq(file1);

    loG("Entering while loop...", Some("info"));

//     let maxkmers = 200;
//     let mut numkmers = 0;
    let mut itrecord : u32 = 0;                                 // We're using it to add the previous indexes!
    let mut theseq   : Vec<u64> = Vec::new();
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        // put_these_nts_into_an_efficient_vector(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);

        // let rl = seqrec.num_bases();
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
            // &itrecord,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                loG(format!("\n\t- COLLISIONS 1 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                loG(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
            }
            // } else {
            //     println!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    loG(format!("\n\t- COLLISIONS 2 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                    loG(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
                }
                // } else {
                //     println!("\n\t\t- NOT COLLISIONS!!!");
                // }
                minmaxdict.entry(hnc).or_insert(hc);

//                 numkmers += 1;
//                 if numkmers >= maxkmers {break};
            }
        }
//         if numkmers >= maxkmers {itrecord += rl as u32;break};
        itrecord += rl as u32;
    }
    loG("Finished getting kmers from first file. Starting with the second...", Some("info"));
//     numkmers = 0;

    let mut reader = open_fastq(file2);

//     println!("MITAD Length of seq. vec.: {}, total length of first files: {}, THING {}, number of recs: {}", theseq.len(), itrecord, (itrecord % 32) as u8, realit);
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
        // put_these_nts_into_an_efficient_vector_rc(&seqrec.seq(), &mut theseq, (itrecord % 32) as u8);
//         println!("{}", seqrec.num_bases());
        // let rl = seqrec.num_bases();
        let rl = seqrec.seq().len();
        let kmer_opt = Kmer::<IntT>::new(
            std::borrow::Cow::Borrowed(seqrec.seq()),
            rl,
            Some(seqrec.qual()),
            k,
            qual.min_qual,
            true,
            // &itrecord,
        );
        if let Some(mut kmer_it) = kmer_opt {
//             numkmers += 1;
            let (hc, hnc, b, km) = kmer_it.get_curr_kmerhash_and_bases_and_kmer();
            outvec.push( (hc, hnc, b) );
            let testkm = outdict.entry(hc).or_insert(km);
            if *testkm != km {
                loG(format!("\n\t- COLLISIONS 3 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                loG(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
            }
            // } else {
            //     println!("\n\t\t- NOT COLLISIONS!!!");
            // }
            minmaxdict.entry(hnc).or_insert(hc);
            while let Some(tmptuple) = kmer_it.get_next_kmer_and_give_us_things() {
                let (hc, hnc, b, km) = tmptuple;
                outvec.push( (hc, hnc, b) );
                let testkm = outdict.entry(hc).or_insert(km);
                if *testkm != km {
                    loG(format!("\n\t- COLLISIONS 4 !!! Hash: {:?}", hc).as_str(), Some("warn"));
                    loG(format!("{:#0258b}\n{:#0258b}", *testkm, km).as_str(), Some("warn"));
                }
                // } else {
                //     println!("\n\t\t- NOT COLLISIONS!!!");
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
//         println!("{:#066b}", tmpu64);
        theseq.push(tmpu64);
    }

    loG("Finished getting kmers from the second file", Some("info"));
    loG(format!("Length of seq. vec.: {}, total length of both files: {}", theseq.len(), itrecord).as_str(), Some("debug"));


    (theseq, outdict, minmaxdict)
}

#[cfg(not(feature = "wasm"))]
fn get_map_with_counts_with_hashes_only(
    invec:      &Vec<(u64, u64, u8)>,
    min_count:  u16,
) -> HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>
{
    let mut outdict = HashMap::with_hasher(BuildHasherDefault::default());


    let mut i = 0;
    let mut c : u16 = 0;
    let mut tmphash = invec[i].0;
    let mut tmpcounter = 0;

    let mut plotvec : Vec<u16> = Vec::new(); // For plotting

    while i < invec.len() {
        if tmphash != invec[i].0 {
            if c >= min_count {
                tmpcounter += 1;
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
//             else {
//                 println!("{}", c);
//             }
//             if i > 500 {
//                 println!("{:?}", outdict);
//                 println!("{} {} {}", tmphash, c, i);
//                 exit(1);
//             }
            tmphash = invec[i].0;
            c = 1;
        } else {
            c = c.saturating_add(1);
        }
        i += 1;
    }


    if c >= min_count {
        tmpcounter += 1;
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

    let root = BitMapBackend::new("./draftrun/histogram.png", (1280, 960)).into_drawing_area();

    let _ = root.fill(&WHITE);

    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .caption("k-mer counts", ("sans-serif", 30.0))
        .build_cartesian_2d((0u32..200u32).into_segmented(), 0u32..200000u32).unwrap();

    chart
        .configure_mesh()
        .disable_x_mesh()
        .bold_line_style(WHITE.mix(0.3))
        .y_desc("Number")
        .x_desc("Count")
        .axis_desc_style(("sans-serif", 15))
        .draw().unwrap();

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.filled())
            .data(plotvec.iter().map(|x: &u16| (*x as u32, 1)).chain(outdict.iter().map(|(_, x)| (x.borrow().counts as u32, 1)))),
    ).unwrap();

    root.present().expect("Unable to write result to file, please make sure 'images' dir exists under current dir");

//     exit(1);

    println!("Good kmers {}", tmpcounter);
    outdict
}


#[cfg(feature = "wasm")]
fn get_map_for_wasm(
    invec:      &Vec<(u64, u64, u8)>,
    min_count:  u16,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u16>)
{
    let mut outdict = HashMap::with_hasher(BuildHasherDefault::default());


    let mut i = 0;
    let mut c : u16 = 0;
    let mut tmphash = invec[i].0;
    let mut tmpcounter = 0;

    let mut plotvec : Vec<u16> = Vec::new(); // For plotting

    while i < invec.len() {
        if tmphash != invec[i].0 {
            if c >= min_count {
                tmpcounter += 1;
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
//             else {
//                 println!("{}", c);
//             }
//             if i > 500 {
//                 println!("{:?}", outdict);
//                 println!("{} {} {}", tmphash, c, i);
//                 exit(1);
//             }
            tmphash = invec[i].0;
            c = 1;
        } else {
            c = c.saturating_add(1);
        }
        i += 1;
    }


    if c >= min_count {
        tmpcounter += 1;
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

    loG(format!("Good kmers {}", tmpcounter).as_str(), Some("debug"));
    (outdict, plotvec)
}


/// Read fastq files, get the reads, get the k-mers, count them, filter them by count, and get some way of recovering the sequence later.
#[cfg(not(feature = "wasm"))]
pub fn preprocessing_gpulike_with_dict_and_seq<IntT>(
    input_files:    &[InputFastx],
    k:              usize,
    qual:           &QualOpts,
    timevec:        &mut Vec<Instant>,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u64>, HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>)
where
    IntT: for<'a> UInt<'a>,
{
    // Build indexes
    log::info!("Starting preprocessing_gpulike with k = {k}");

    if any_fastq(input_files) {
        log::info!("FASTQ files filtered with: {qual}");
    } else {
        log::info!("All input files FASTA (no error filtering)");
    }

    let total_size  = input_files.len();
    if total_size > 1 {panic!("Not expecting more than a pair of pair-end reads right now!");};


    // First, we want to fill our mega-vector with all k-mers from both paired-end reads
    log::info!("Filling vector");

    let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
    let (theseq, thedict, maxmindict) = get_kmers_from_both_files_and_the_dict_and_the_seq::<IntT>(&input_files[0].1,
        input_files[0].2.as_ref().expect("No paired reads!"),
        k,
        qual,
        &mut tmpvec);

    timevec.push(Instant::now());
    log::info!("Kmers extracted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());


    // Then, we want to sort it according to the hash
    println!("Number of kmers BEFORE cleaning: {:?}", tmpvec.len());
    log::info!("Sorting vector");
    tmpvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    timevec.push(Instant::now());
    log::info!("Kmers sorted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());


    // Then, do a counting of everything and save the results in a dictionary and return it
    timevec.push(Instant::now());
    log::info!("Counting k-mers");
    let themap = get_map_with_counts_with_hashes_only(&tmpvec, qual.min_count);
    drop(tmpvec);

    timevec.push(Instant::now());
    log::info!("Kmers counted in {} s", timevec.last().unwrap().duration_since(*timevec.get(timevec.len().wrapping_sub(2)).unwrap()).as_secs());

    (themap, theseq, thedict, maxmindict)
}



#[cfg(feature = "wasm")]
pub fn preprocessing_for_wasm<IntT>(
    file1:    &mut WebSysFile,
    file2:    &mut WebSysFile,
    k:         usize,
    qual:     &QualOpts,
) -> (HashMap::<u64, RefCell<HashInfoSimple>, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u64>, Option<HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>>, HashMap::<u64, u64, BuildHasherDefault<NoHashHasher<u64>>>, Vec<u16>)
where
    IntT: for<'a> UInt<'a>,
{
    // Build indexes
    loG("Starting preprocessing with k = {k}", Some("info"));


    // First, we want to fill our mega-vector with all k-mers from both paired-end reads
    loG("Filling vector", Some("info"));

    let mut tmpvec : Vec<(u64, u64, u8)> = Vec::new();
    let (theseq, thedict, maxmindict) = get_kmers_from_both_files_wasm::<IntT>(file1,
        file2,
        k,
        qual,
        &mut tmpvec);

    loG("Kmers extracted", Some("info"));


    // Then, we want to sort it according to the hash
    // println!("Number of kmers BEFORE cleaning: {:?}", tmpvec.len());
    loG("Sorting vector", Some("info"));
    tmpvec.par_sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    loG("Kmers sorted.", Some("info"));


    // Then, do a counting of everything and save the results in a dictionary and return it
    loG("Counting k-mers", Some("info"));
    let (themap, histovec) = get_map_for_wasm(&tmpvec, qual.min_count);
    drop(tmpvec);

    loG("Kmers counted.", Some("info"));

    (themap, theseq, Some(thedict), maxmindict, histovec)
}
