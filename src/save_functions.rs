//! Some docs should be here

use core::panic;
use nohash_hasher::NoHashHasher;
use std::{collections::HashMap, hash::BuildHasherDefault};

use super::io_utils::*;
// use std::process::exit;

use crate::bit_encoding::UInt;
use crate::graph_works::Contigs;



/// Writes the contig sequences and hopefully their average counts/coverage in the future
pub fn write_sequences_and_coverages<IntT>(invec : &mut Contigs,
                                           inmap : &    HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>,
                                           k     :      usize,)
where
    IntT: for<'a> UInt<'a>, {
    // TODO: implement coverages somehow...
    invec.contig_sequences = Some(Vec::with_capacity(invec.serialized_contigs.len()));
    let mut counter = 0;
    for ipc in 0..invec.serialized_contigs.len() {
//         println!("\nIteration");
//         let mut outseq : VecDeque<u8>  = VecDeque::new();
        let mut outseq : Vec<u8>  = Vec::new();

//         println!("Initial index: {}", initind);
        // First of all, we decode and set the cov. for the first k nucleotides
        if invec.serialized_contigs[ipc].len() > 1 {panic!("MORE THAN ONE ENTRY!!")};
        let initkmer    = inmap.get(&invec.serialized_contigs[ipc][0].abs_ind[0]).unwrap();
        let nbitstomove = initkmer.n_bits() as usize - 2 * (k - 1);
//         println!("nbits: {nbitstomove}");

        let mut prevkmer = *initkmer;
        if invec.serialized_contigs[ipc][0].abs_ind.len() > 1 {
            let tmpkmermoved    = *inmap.get(&invec.serialized_contigs[ipc][0].abs_ind[1]).unwrap() >> 2;
            let tmpkmerrevmoved = inmap.get(&invec.serialized_contigs[ipc][0].abs_ind[1]).unwrap().rev_comp(k) >> 2;
            let prevkmermoved   = (prevkmer << nbitstomove) >> nbitstomove;
//             println!("Prev.              {:#066b}", prevkmer            );
//             println!("Prev. (rev.-comp.) {:#066b}", prevkmer.rev_comp(k));
//             println!("Post.              {:#066b}", tmpkmermoved);
//             println!("Post. (rev.-comp.) {:#066b}", tmpkmerrevmoved);
//             println!("Prev.              {:#066b}", prevkmermoved);

            if tmpkmermoved != prevkmermoved && tmpkmerrevmoved != prevkmermoved {
                // We need to add the first nucleotides from the rev. comp.
                prevkmer = prevkmer.rev_comp(k);
//                 println!("CHANGED");
            }

            if tmpkmermoved == prevkmermoved && tmpkmerrevmoved == prevkmermoved {
                panic!("HOLI");
            }
        }

        // for inc in 0..k {
        //     outseq.push( prevkmer.get_one_nucleotide(k - 1 - inc));
        // }

        outseq.push( prevkmer.get_one_nucleotide(0));

        // And now, we start the hard work with the remaining nucleotides
        let mut currkmer : IntT;
        let thelen = invec.serialized_contigs[ipc][0].abs_ind.len();
        for i in 1..thelen {
//             println!("Entry {}/{}", i + 1, thelen);
            currkmer = *inmap.get(&invec.serialized_contigs[ipc][0].abs_ind[i]).unwrap();

            if ((prevkmer << nbitstomove) >> nbitstomove) == (currkmer >> 2) && ((prevkmer << nbitstomove) >> nbitstomove) == (currkmer.rev_comp(k) >> 2) {
                panic!("HOLI");
            }

//             if ((prevkmer.rev_comp(k) << nbitstomove) >> nbitstomove) == (currkmer >> 2) || ((prevkmer.rev_comp(k) << nbitstomove) >> nbitstomove) == currkmer.rev_comp(k) {
//                 panic!("TEST2");
//             }

            if ((prevkmer << nbitstomove) >> nbitstomove) != (currkmer >> 2) {
//                 println!("CAMBIANDO!");
                currkmer = currkmer.rev_comp(k);
                if ((prevkmer << (nbitstomove)) >> nbitstomove) != (currkmer >> 2) {
//                     println!("BAD THING");
//                     println!("Prev.:              {:#066b}", (prevkmer << (nbitstomove)) >> nbitstomove);
//                     println!("Prev. (rev.-comp.): {:#066b}", (prevkmer.rev_comp(k) << (nbitstomove)) >> nbitstomove);
//                     println!("Post:               {:#066b}", currkmer.rev_comp(k) >> 2);
//                     println!("Post. (rev.-comp.): {:#066b}", currkmer >> 2);
//                     println!("Prev. hash: {}",   invec.serialized_contigs[ipc][0].abs_ind[i - 1]);
//                     println!("Post. hash: {}\n", invec.serialized_contigs[ipc][0].abs_ind[i]);
                    counter += 1;
                }
            }

            outseq.push(currkmer.get_one_nucleotide(0));
            prevkmer = currkmer;

        }

        if outseq.len() > k {
            for _ in 0..k {
                outseq.remove(outseq.len() - 1);
            }
            if outseq.len() > 100 {
                invec.contig_sequences.as_mut().unwrap().push(outseq);
            }
        }

    }
    println!("\nNUMBER OF BAD THINGS: {}\n", counter);

}


/// Stores all the contigs as a fasta file
#[cfg(not(feature = "wasm"))]
pub fn save_as_fasta<IntT>(ingraph: &mut Contigs,
                                            inmap:   &    HashMap::<u64, IntT, BuildHasherDefault<NoHashHasher<u64>>>,
                                            k:        usize,
                                            outfile: &String)
where
    IntT: for<'a> UInt<'a>, {
    // First, we write the sequences and the coverages
    write_sequences_and_coverages(ingraph, inmap, k);

    println!("Starting to save");
    println!("{}", outfile);
    // Now, we just write all the contigs. We get our writing buffer with this:
    let mut wbuf = set_ostream(&Some(outfile.clone()));
    // And simply, contig per contig, we write the file

    println!("\tLen.\tMean\tSD\tMedian");
    ingraph.write_fasta(&mut wbuf);
}
