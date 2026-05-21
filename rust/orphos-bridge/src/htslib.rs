pub mod bgzf;
pub mod faidx;
pub mod tabix;

use std::io::Cursor;

pub struct FaidxResult {
    pub fai: Vec<u8>,
    pub gzi: Vec<u8>,
}

impl FaidxResult {
    pub fn take_fai(self) -> Vec<u8> {
        self.fai
    }
    pub fn take_gzi(self) -> Vec<u8> {
        self.gzi
    }
}

/// Compress `data` into a BGZF stream and return the bytes.
pub fn bgzf_compress(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    bgzf::bgzf_compress(Cursor::new(data), &mut out).expect("bgzf compress failed");
    out
}

/// Index a BGZF-compressed FASTA; returns FAI and GZI bytes.
pub fn faidx_index_fasta(bgzf_data: &[u8]) -> FaidxResult {
    let mut fai = Vec::new();
    let mut gzi = Vec::new();
    faidx::faidx_index_fasta(Cursor::new(bgzf_data), &mut fai, &mut gzi)
        .expect("faidx index failed");
    FaidxResult { fai, gzi }
}

/// Build a CSI index for a BGZF-compressed GFF3; returns the BGZF-compressed CSI bytes.
pub fn csi_index_gff(bgzf_data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    tabix::csi_index_gff(Cursor::new(bgzf_data), &mut out).expect("csi index failed");
    out
}
