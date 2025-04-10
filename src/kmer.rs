//! Some docs here, but I think the file name is already self-explicative enough
//!

use std::borrow::Cow;

use super::bit_encoding::*;
use super::nthash::NtHashIterator;

// use std::process::exit;


/// Struct to generate all k-mers from an input sequence
///
/// Holds reference to input sequence, current encoded k-mer, and other
/// information (k-mer size, masks etc.)
#[derive(Debug)]
pub struct Kmer<'a, IntT> {
    /// K-mer size
    k: usize,
    /// Reference to input sequence
    seq: Cow<'a, [u8]>,
    /// Size of seq
    seq_len: usize,
    /// Reference to sequence quality scores
    qual: Option<&'a [u8]>,
    /// Minimum quality score to allow in a k-mer
    min_qual: u8,
    /// Current index in input sequence
    index: usize,
    /// Current k-mer
    kmer: IntT,
    /// Hash generator for reads
    hash_gen: NtHashIterator,
    /// rc or Not
    rc: bool,
    // /// the index of the record
    // rec_ind: &'a u32,
}


impl<'a, IntT: for<'b> UInt<'b>> Kmer<'a, IntT> {
    /// Quality score is at least minimum.
    #[inline(always)]
    fn valid_qual(idx: usize, qual: Option<&'a [u8]>, min_qual: u8) -> bool {
        match qual {
            Some(qual_seq) => (qual_seq[idx] - 33) > min_qual, // ASCII encoding starts from b'!' = 33
            None => true,
        }
    }

    /// Build a new k-mer at the given index.
    ///
    #[allow(clippy::too_many_arguments)]
    fn build(
        seq: &[u8],
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        idx: &mut usize,
        min_qual: u8,
        rc: bool,
        // rec_ind: &'a u32,
    ) -> Option<(IntT, NtHashIterator)> {
        // log::info!("Building kmer");
        if *idx + k >= seq_len {
            return None;
        }
        let mut kmer = IntT::zero_init();
        let mut i = 0;
        while i < k {
            if valid_base(seq[i + *idx]) && (Self::valid_qual(i + *idx, qual, min_qual)) { // Checks for N or n
                let next_base = encode_base(seq[i + *idx]);
                kmer <<= 2;
                // println!("Adding {:#0258b}", IntT::from_encoded_base(next_base));
                kmer |= IntT::from_encoded_base(next_base);
                i += 1;
            } else {
                // Start again, skipping over N
                // println!("Bad base found, restarting!");
                *idx += i + 1;
                if *idx + k >= seq_len {
                    return None;
                }
                kmer = IntT::zero_init();
                i    = 0;
            }
        }
        // println!("Final kmer:\t {:#0258b}", kmer);
        // println!("{:#064b} {:#064b} {:#064b}", ((kmer).as_u8() & 3), ((kmer >> (k - 1)*2) << 2).as_u8(), ((kmer).as_u8() & 3) | ((kmer >> (k - 1)*2) << 2).as_u8());
        let hash_gen = NtHashIterator::new(&seq[*idx..(*idx + k)], k, rc);
        *idx += k - 1;
//         log::info!("kmer built!");
        Some((kmer, hash_gen))
    }

    /// Move forward to the next valid k-mer
    ///
    /// Usually the next base, but if an N skips over it.
    /// If end of sequence encountered then returns `false`.
    fn roll_fwd(&mut self) -> bool {
//         log::info!("Rolling kmer");
        let mut success = false;
        self.index += 1;
        if self.index >= self.seq_len {
            return success;
        }
        let base = self.seq[self.index];
        if !valid_base(base) || !Self::valid_qual(self.index, self.qual, self.min_qual) {
            // log::info!("Next base is NOT valid, we get a new kmer");
            let new_kmer = Self::build(
                &self.seq,
                self.seq_len,
                self.qual,
                self.k,
                &mut self.index,
                self.min_qual,
                self.rc,
                // self.rec_ind,
            );
            if let Some(kmer_tuple) = new_kmer {
                (self.kmer, self.hash_gen) = kmer_tuple;
                success = true;
            }
        } else {
            // log::info!("Next base is valid, we're rolling, baby!");
            let new_base = encode_base(base);
//             println!("");
            // println!("cur kmer:\t {:#0258b}", self.kmer);

            let old_base = (self.kmer >> ((self.k - 1) * 2)).as_u8() & 3;
            // println!("new base:\t {:#0258b}", new_base);
            // println!("old base:\t {:#0258b}", old_base);
            self.hash_gen.roll_fwd(old_base, new_base);

            // Update the k-mer
            // println!("previous kmer:\t {:#0258b}", self.kmer);
            let cleanbits = self.kmer.n_bits() as usize - (self.k - 1) * 2;
            self.kmer = (((self.kmer << cleanbits) >> cleanbits) << 2) | (IntT::from_encoded_base(new_base)) ;
            // println!("new kmer:\t {:#0258b}\n", self.kmer);
            success = true;
        }
        // log::info!("Finishing rolling");
        success
    }


    /// Create a [`Kmer`] iterator given reference to sequence input.
    ///
    /// Sequence, length and quality come from [`needletail`].
    ///
    /// Returns [`None`] if no valid k-mers found in input (e.g. too short,
    /// no sequence, too many Ns).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        seq: Cow<'a, [u8]>,
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        min_qual: u8,
        rc: bool,
        // rec_ind: &'a u32,
    ) -> Option<Self> {
        let mut index = 0;
        let first_kmer = Self::build(
            &seq,
            seq_len,
            qual,
            k,
            &mut index,
            min_qual,
            rc,
            // rec_ind,
        );
        if let Some((kmer, hash_gen)) = first_kmer {
            let kmerobj = Self {
                k,
                seq,
                seq_len,
                qual,
                min_qual,
                index,
                kmer,
                hash_gen,
                rc,
                // rec_ind,
            };
            Some(kmerobj)
        } else {
            None
        }
    }

    /// Get the current k-mer hash
    ///
    /// Returns the canonical k-mer hash, the non-canonical one, the first and last bases, and the k-mer sequence
    pub fn get_curr_kmerhash_and_bases_and_kmer(&self) -> (u64, u64, u8, IntT) {
        // OLD bases extraction: we need now to see if we are
        let (canhash, notcanhash, isittherevcomp) = self.hash_gen.curr_hash_and_whether_it_is_the_inverse();
        let thebases : u8; // Last four bits correspond to the canonical kmer bases, first four bits to the alternative kmer bases
//         if isittherevcomp {
//             thebases = (rc_base((self.kmer >> (self.k - 1)*2).as_u8()) & 3)
//                         | (rc_base((self.kmer).as_u8() & 3) << 2)
//                         | ((self.kmer >> (self.k - 1)*2) << 6).as_u8()
//                         | (((self.kmer).as_u8() & 3) << 4);
//         } else {
//             thebases =  (((self.kmer >> (self.k - 1)*2).as_u8() & 3) << 2)
//                         | ((self.kmer).as_u8() & 3)
//                         | (rc_base((self.kmer).as_u8() & 3) << 6)
//                         | (rc_base((self.kmer >> (self.k - 1)*2).as_u8() & 3) << 4);
//         }
//             thebases =  (((self.kmer >> (self.k - 1)*2).as_u8() & 3) << 2)
//                         | ((self.kmer).as_u8() & 3);

        let thekmer : IntT;
//         println!("{} {} {}", canhash, notcanhash, isittherevcomp);
//         println!("{:#0194b}\n{:#0194b}\n", self.kmer, self.kmer.rev_comp(self.k));
        if isittherevcomp {
            thebases = (rc_base((self.kmer >> (self.k - 1)*2).as_u8()) & 3) | (rc_base((self.kmer).as_u8() & 3) << 2);
            thekmer  = self.kmer.rev_comp(self.k);
        } else {
            thebases = ((self.kmer >> (self.k - 1)*2) << 2).as_u8() | ((self.kmer).as_u8() & 3);
            thekmer  = self.kmer;
        }
//         println!("{:#010b}\n", thebases);

//         println!("{:#066b} {}", self.kmer.to_u64(), self.get_hash());
//         println!("{:#066b}", thebases);

        (canhash, notcanhash, thebases, thekmer)
    }

    /// Get a `u64` hash of the current k-mer using [`NtHashIterator`]
    ///
    /// # Panics
    ///
    /// If called after creating with `is_reads = false`
    pub fn get_hash(&self) -> u64 {
        self.hash_gen.curr_hash()
    }

    // /// Get the next k-mer in the sequence
    // ///
    // pub fn get_next_kmer(&mut self) -> Option<(u64, u8, u64)> {
    //     let next = self.roll_fwd();
    //     match next {
    //         true => Some(self.get_curr_kmerhash_and_bases()),
    //         false => None,
    //     }
    // }

    /// Get the next k-mer in the sequence and provide other information
    pub fn get_next_kmer_and_give_us_things(&mut self) -> Option<(u64, u64, u8, IntT)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmerhash_and_bases_and_kmer()),
            false => None,
        }
    }

    // pub fn create_absolute_index_for_illumina(&self) -> u64 {
    // //     println!("{} {} {}", *seq_ind, *rec_ind, *rec_ind as u64 + *seq_ind as u64);
    //     if self.rc {
    //         return *self.rec_ind as u64 + ((self.seq_len - self.index - 1 + self.k - 1) as u64);
    //     } else {
    //         return *self.rec_ind as u64 + self.index as u64;
    //     }
    // }
}

// fn create_absolute_index(rc: &bool,
//     seq_ind: &usize,
//     rec_ind: u32) -> u64 {
//
//     let mut outind : u64 = (rec_ind as u64) << 1;
//     if *rc {
//         outind |= 1u64;
//     }
//
//     outind <<= 32;
//
//     outind | ((*seq_ind as u64) & 33u64)
// }


