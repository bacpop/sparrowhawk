//! ntHash is a hash function optimised for a DNA alphabet `{A, C, T, G}`.
//!
//! It works particularly well as a rolling hash, e.g. for k-mers in an input
//! sequence.
//!
//! This implementation based on ntHash 2

use super::bit_encoding::{encode_base, rc_base};

/// Values to look up and multiply for direct k-mers
pub const HASH_LOOKUP: [u64; 4] = [
    0x3c8b_fbb3_95c6_0474,          // A
    0x3193_c185_62a0_2b4c,          // C
    0x2955_49f5_4be2_4456,          // T
    0x2032_3ed0_8257_2324,          // G
];

/// Values to look up and multiply for reverse-complement k-mers
pub const RC_HASH_LOOKUP: [u64; 4] = [
    0x2955_49f5_4be2_4456,
    0x2032_3ed0_8257_2324,
    0x3c8b_fbb3_95c6_0474,
    0x3193_c185_62a0_2b4c,
];

/// This is A33R, C33R, T33R, G33R so can be indexed by nucleotide * 33
pub const MS_TAB_33R: [u64; 132] = [
    0x195c60474,
    0x12b8c08e9,
    0x571811d3,
    0xae3023a6,
    0x15c60474c,
    0xb8c08e99,
    0x171811d32,
    0xe3023a65,
    0x1c60474ca,
    0x18c08e995,
    0x11811d32b,
    0x3023a657,
    0x60474cae,
    0xc08e995c,
    0x1811d32b8,
    0x1023a6571,
    0x474cae3,
    0x8e995c6,
    0x11d32b8c,
    0x23a65718,
    0x474cae30,
    0x8e995c60,
    0x11d32b8c0,
    0x3a657181,
    0x74cae302,
    0xe995c604,
    0x1d32b8c08,
    0x1a6571811,
    0x14cae3023,
    0x995c6047,
    0x132b8c08e,
    0x6571811d,
    0xcae3023a,
    0x162a02b4c,
    0xc5405699,
    0x18a80ad32,
    0x115015a65,
    0x2a02b4cb,
    0x54056996,
    0xa80ad32c,
    0x15015a658,
    0xa02b4cb1,
    0x140569962,
    0x80ad32c5,
    0x1015a658a,
    0x2b4cb15,
    0x569962a,
    0xad32c54,
    0x15a658a8,
    0x2b4cb150,
    0x569962a0,
    0xad32c540,
    0x15a658a80,
    0xb4cb1501,
    0x169962a02,
    0xd32c5405,
    0x1a658a80a,
    0x14cb15015,
    0x9962a02b,
    0x132c54056,
    0x658a80ad,
    0xcb15015a,
    0x1962a02b4,
    0x12c540569,
    0x58a80ad3,
    0xb15015a6,
    0x14be24456,
    0x97c488ad,
    0x12f89115a,
    0x5f1222b5,
    0xbe24456a,
    0x17c488ad4,
    0xf89115a9,
    0x1f1222b52,
    0x1e24456a5,
    0x1c488ad4b,
    0x189115a97,
    0x11222b52f,
    0x24456a5f,
    0x488ad4be,
    0x9115a97c,
    0x1222b52f8,
    0x4456a5f1,
    0x88ad4be2,
    0x1115a97c4,
    0x22b52f89,
    0x456a5f12,
    0x8ad4be24,
    0x115a97c48,
    0x2b52f891,
    0x56a5f122,
    0xad4be244,
    0x15a97c488,
    0xb52f8911,
    0x16a5f1222,
    0xd4be2445,
    0x1a97c488a,
    0x152f89115,
    0xa5f1222b,
    0x82572324,
    0x104ae4648,
    0x95c8c91,
    0x12b91922,
    0x25723244,
    0x4ae46488,
    0x95c8c910,
    0x12b919220,
    0x57232441,
    0xae464882,
    0x15c8c9104,
    0xb9192209,
    0x172324412,
    0xe4648825,
    0x1c8c9104a,
    0x191922095,
    0x12324412b,
    0x46488257,
    0x8c9104ae,
    0x11922095c,
    0x324412b9,
    0x64882572,
    0xc9104ae4,
    0x1922095c8,
    0x124412b91,
    0x48825723,
    0x9104ae46,
    0x122095c8c,
    0x4412b919,
    0x88257232,
    0x1104ae464,
    0x2095c8c9,
    0x412b9192,
];

/// This is A31R, C31R, T31R, G31R so can be indexed by nucleotide * 31
pub const MS_TAB_31L: [u64; 124] = [
    0x3c8bfbb200000000,
    0x7917f76400000000,
    0xf22feec800000000,
    0xe45fdd9200000000,
    0xc8bfbb2600000000,
    0x917f764e00000000,
    0x22feec9e00000000,
    0x45fdd93c00000000,
    0x8bfbb27800000000,
    0x17f764f200000000,
    0x2feec9e400000000,
    0x5fdd93c800000000,
    0xbfbb279000000000,
    0x7f764f2200000000,
    0xfeec9e4400000000,
    0xfdd93c8a00000000,
    0xfbb2791600000000,
    0xf764f22e00000000,
    0xeec9e45e00000000,
    0xdd93c8be00000000,
    0xbb27917e00000000,
    0x764f22fe00000000,
    0xec9e45fc00000000,
    0xd93c8bfa00000000,
    0xb27917f600000000,
    0x64f22fee00000000,
    0xc9e45fdc00000000,
    0x93c8bfba00000000,
    0x27917f7600000000,
    0x4f22feec00000000,
    0x9e45fdd800000000,
    0x3193c18400000000,
    0x6327830800000000,
    0xc64f061000000000,
    0x8c9e0c2200000000,
    0x193c184600000000,
    0x3278308c00000000,
    0x64f0611800000000,
    0xc9e0c23000000000,
    0x93c1846200000000,
    0x278308c600000000,
    0x4f06118c00000000,
    0x9e0c231800000000,
    0x3c18463200000000,
    0x78308c6400000000,
    0xf06118c800000000,
    0xe0c2319200000000,
    0xc184632600000000,
    0x8308c64e00000000,
    0x6118c9e00000000,
    0xc23193c00000000,
    0x1846327800000000,
    0x308c64f000000000,
    0x6118c9e000000000,
    0xc23193c000000000,
    0x8463278200000000,
    0x8c64f0600000000,
    0x118c9e0c00000000,
    0x23193c1800000000,
    0x4632783000000000,
    0x8c64f06000000000,
    0x18c9e0c200000000,
    0x295549f400000000,
    0x52aa93e800000000,
    0xa55527d000000000,
    0x4aaa4fa200000000,
    0x95549f4400000000,
    0x2aa93e8a00000000,
    0x55527d1400000000,
    0xaaa4fa2800000000,
    0x5549f45200000000,
    0xaa93e8a400000000,
    0x5527d14a00000000,
    0xaa4fa29400000000,
    0x549f452a00000000,
    0xa93e8a5400000000,
    0x527d14aa00000000,
    0xa4fa295400000000,
    0x49f452aa00000000,
    0x93e8a55400000000,
    0x27d14aaa00000000,
    0x4fa2955400000000,
    0x9f452aa800000000,
    0x3e8a555200000000,
    0x7d14aaa400000000,
    0xfa29554800000000,
    0xf452aa9200000000,
    0xe8a5552600000000,
    0xd14aaa4e00000000,
    0xa295549e00000000,
    0x452aa93e00000000,
    0x8a55527c00000000,
    0x14aaa4fa00000000,
    0x20323ed000000000,
    0x40647da000000000,
    0x80c8fb4000000000,
    0x191f68200000000,
    0x323ed0400000000,
    0x647da0800000000,
    0xc8fb41000000000,
    0x191f682000000000,
    0x323ed04000000000,
    0x647da08000000000,
    0xc8fb410000000000,
    0x91f6820200000000,
    0x23ed040600000000,
    0x47da080c00000000,
    0x8fb4101800000000,
    0x1f68203200000000,
    0x3ed0406400000000,
    0x7da080c800000000,
    0xfb41019000000000,
    0xf682032200000000,
    0xed04064600000000,
    0xda080c8e00000000,
    0xb410191e00000000,
    0x6820323e00000000,
    0xd040647c00000000,
    0xa080c8fa00000000,
    0x410191f600000000,
    0x820323ec00000000,
    0x40647da00000000,
    0x80c8fb400000000,
    0x10191f6800000000,
];


/// Stores forward and (optionally) reverse complement hashes of k-mers in a nucleotide sequence
#[derive(Debug)]
pub struct NtHashIterator {
    k: usize,
    fh: u64,
    rh: Option<u64>,
}


/// Function that exchanges the bits in the 0 and 33 positions
#[inline(always)]
pub fn swapbits033(v: u64) -> u64 {
    let x = (v ^ (v >> 33)) & 1;
    v ^ (x | (x << 33))
}


/// Function that exchanges the bits in the 32 and 63 positions
#[inline(always)]
pub fn swapbits3263(v: u64) -> u64 {
    let x = ((v >> 32) ^ (v >> 63)) & 1;
    v ^ ((x << 32) | (x << 63))
}


impl NtHashIterator {
    /// Creates a new iterator over a sequence with a given k-mer size
    pub fn new(seq: &[u8], k: usize, rc: bool) -> NtHashIterator {
        let mut fh : u64 = 0;
        // for (_, v) in seq[0..k].iter().enumerate() {
        for v in seq[0..k].iter() {
            fh = fh.rotate_left(1_u32); // This, WITH swapbits033 (next line), is the "srol" function
            fh = swapbits033(fh);
            fh ^= HASH_LOOKUP[encode_base(*v) as usize];
        }

        let rh = if rc {
            let mut h : u64 = 0;
            // for (_, v) in seq[0..k].iter().rev().enumerate() {
            for v in seq[0..k].iter().rev() {
                h =  h.rotate_left(1_u32);
                h =  swapbits033(h);
                h ^= RC_HASH_LOOKUP[encode_base(*v) as usize];
            }
            Some(h)
        } else {
            None
        };

        Self { k, fh, rh }
    }

    /// Move to the next k-mer by adding a new base, removing a base from the end, efficiently updating the hash.
    pub fn roll_fwd(&mut self, old_base: u8, new_base: u8) {
        self.fh = self.fh.rotate_left(1);
        self.fh = swapbits033(self.fh);
        self.fh ^= HASH_LOOKUP[new_base as usize];
        self.fh ^=   MS_TAB_31L[(old_base as usize * 31) + (self.k % 31)]
                   | MS_TAB_33R[(old_base as usize) * 33 + (self.k % 33)];



        if let Some(rev) = self.rh {
            let mut h = rev
                ^ (MS_TAB_31L[(rc_base(new_base) as usize * 31) + (self.k % 31)]
                 | MS_TAB_33R[(rc_base(new_base) as usize) * 33 + (self.k % 33)]);
            h ^= RC_HASH_LOOKUP[old_base as usize];
            h = h.rotate_right(1_u32);
            h = swapbits3263(h);
            self.rh = Some(h);
        };
    }

    /// Retrieve the current hash (minimum of forward and reverse complement hashes)
    pub fn curr_hash(&self) -> u64 {
        if let Some(rev) = self.rh {
            u64::min(self.fh, rev)
        } else {
            self.fh
        }
    }

    /// Retrieve the current hash (minimum of forward and reverse complement hashes)
    pub fn curr_hash_and_whether_it_is_the_inverse(&self) -> (u64, u64, bool) {
        let rev = self.rh.unwrap();
        (u64::min(self.fh, rev), u64::max(self.fh, rev), if rev < self.fh {true} else {false})
    }
}
