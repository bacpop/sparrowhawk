//! A Bloom filter for counting k-mers from reads
//!
//! This is to support a minimum count from FASTQ files, for basic error
//! correction, while keeping memory use low and predictable. Use the
//! [`KmerFilter`] struct.

use std::borrow::BorrowMut;
use std::cmp::Ordering;

use std::{collections::HashMap, hash::BuildHasherDefault};
use nohash_hasher::NoHashHasher;


/// Default bloom filter width (expected number of k-mers)
///
/// 2^27 =~ 130M
const BLOOM_WIDTH: usize = 1 << 27;
/// Number of bits to use in each bloom block (~1% FPR)
const BITS_PER_ENTRY: usize = 12;

/// A filter which counts input k-mers, returns whether they have passed a count threshold.
///
/// Uses a blocked bloom filter as a first pass to remove singletons.
/// Code for blocked bloom filter based on:
/// <https://github.com/lemire/Code-used-on-Daniel-Lemire-s-blog/blob/master/2021/10/02/wordbasedbloom.cpp>
///
/// This has the advantage of using less memory than a larger countmin filter,
/// being a bit faster (bloom is ~3x faster than countmin, but having count also
/// allows entry to dictionary to be only checked once for each passing k-mer)
///
/// Once passed through the bloom filter, a HashMap is used for counts >=2.
/// This filter therefore has no false-negatives and negligible false-negatives
#[derive(Debug, Clone, Default)]
pub struct KmerFilter {
    /// Size of the bloom filter
    buf_size: u64,
    /// Buffer for the bloom filter
    buffer: Vec<u64>,
    /// Table of counts
    counts: HashMap::<u64, (u16, u64, u8), BuildHasherDefault<NoHashHasher<u64>>>,
    /// Minimum count to pass filter
    min_count: u16,
}

impl KmerFilter {
    /// Cheap modulo
    /// https://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    #[inline(always)]
    fn reduce(key: u64, range: u64) -> u64 {
        (((key as u128) * (range as u128)) >> 64) as u64
    }

    /// Like splitmix64 but simpler and faster
    #[inline(always)]
    fn cheap_mix(key: u64) -> u64 {
        (key ^ (key >> 31)).wrapping_mul(0x85D0_59AA_3331_21CF)
    }

    /// Extract five bits from word as fingerprint
    #[inline(always)]
    fn fingerprint(key: u64) -> u64 {
        1 << (key & 63)
            | 1 << ((key >> 6) & 63)
            | 1 << ((key >> 12) & 63)
            | 1 << ((key >> 18) & 63)
            | 1 << ((key >> 24) & 63)
    }

    /// Generate a location in the buffer from the hash
    #[inline(always)]
    fn location(key: u64, range: u64) -> usize {
        Self::reduce(Self::cheap_mix(key), range) as usize
    }

    /// Check if in the bloom filter, add if not. Returns whether passed filter
    fn bloom_add_and_check(&mut self, key: u64) -> bool {
        let f_print = Self::fingerprint(key);
        let buf_val = self.buffer[Self::location(key, self.buf_size)].borrow_mut();
        if *buf_val & f_print == f_print {
            true
        } else {
            *buf_val |= f_print;
            false
        }
    }

    /// Creates a new filter with given threshold
    ///
    /// Note:
    /// - Maximum count is [`u16::MAX`] i.e. 65535
    /// - Must call [`KmerFilter::init()`] before using.
    pub fn new(min_count: u16) -> Self {
        let buf_size =
            f64::round(BLOOM_WIDTH as f64 * (BITS_PER_ENTRY as f64 / 8.0) / (u64::BITS as f64))
                as u64;
        Self {
            buf_size,
            buffer: Vec::new(),
            counts: HashMap::with_hasher(BuildHasherDefault::default()),
            min_count,
        }
    }

    /// Initialises table so it is ready for use.
    ///
    /// Allocates memory for bloom filter (which can be avoided with FASTA input)
    pub fn init(&mut self) {
        if self.buffer.is_empty() {
            self.buffer.resize(self.buf_size as usize, 0);
        }
    }

    /// Add an observation of a k-mer and middle base to the filter, and return if it passed
    /// minimum count filtering criterion.
    pub fn filter(&mut self, kmer_hash : u64, kmer_nc_hash : u64, bases : u8) -> Ordering {
        // This is possible because of the k-mer size restriction, the top two
        // bit are always zero
        match self.min_count {
            // No filtering
            0 | 1 => Ordering::Equal,
            // Just the bloom filter
            2 => {
                if self.bloom_add_and_check(kmer_hash) {
                    Ordering::Equal
                } else {
                    Ordering::Less
                }
            }
            // Bloom filter then hash table
            _ => {
                if self.bloom_add_and_check(kmer_hash) {
                    let mut count: u16 = 2;
                    self.counts
                        .entry(kmer_hash)
                        .and_modify(|tuple| {
                            count = tuple.0.saturating_add(1);
                            tuple.0 = count;
                        })
                        .or_insert((count, kmer_nc_hash, bases));
                    self.min_count.cmp(&count)
                } else {
                    Ordering::Less
                }
            }
        }
    }

    /// Get method to retrieve the count map
    pub fn get_counts_map(&mut self) -> &mut HashMap::<u64, (u16, u64, u8), BuildHasherDefault<NoHashHasher<u64>>> {
        return &mut self.counts;
    }
}
