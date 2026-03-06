use std::collections::HashMap;
use std::io::{self, Read, Write};
use super::bgzf::{BgzfReader, BgzfWriter};

// -----------------------------------------------------------------
// CSI format constants (tabix -C -p gff, htslib default)
// -----------------------------------------------------------------
const MIN_SHIFT: u32 = 14;
/// Number of index levels used by `tabix -C` for GFF (= 8),
/// Covers coordinates up to 2^(14+3*8) = 2^38 ≈ 274 GB.
const N_LVLS: u32 = 8;
/// Pseudo-bin for per-sequence metadata (N_BINS + 1).
const META_BIN: u32 = 19_173_962;

/// Minimum compressed-byte span for a bin to be kept at its level rather than
/// merged into its parent (= HTS_MIN_MARKER_DIST = 0x10000 = one BGZF block).
const HTS_MIN_MARKER_DIST: u64 = 0x10000;

// -----------------------------------------------------------------
// Binning helpers
// -----------------------------------------------------------------

/// First bin number at level l.
fn hts_bin_first(l: u32) -> u32 {
    ((1u32 << (3 * l)) - 1) / 7
}

/// Parent bin of b (htslib: (b-1) >> 3).
fn hts_bin_parent(b: u32) -> u32 {
    (b - 1) >> 3
}

/// Level of bin b: count parent steps from b until reaching 0.
fn hts_bin_level(b: u32) -> u32 {
    let mut b = b;
    let mut level = 0;
    while b > 0 {
        b = (b - 1) >> 3;
        level += 1;
    }
    level
}

/// Bottom linear-index slot covered by bin b (with N_LVLS levels).
fn hts_bin_bot(b: u32) -> u64 {
    let level = hts_bin_level(b);
    let offset = b - hts_bin_first(level);
    (offset as u64) << ((N_LVLS - level) * 3)
}

/// Compute loff for a bin: lidx[hts_bin_bot(bin)], falling back to the last
/// non-zero lidx entry (mirrors htslib update_loff).
fn compute_loff(bin: u32, lidx: &[u64]) -> u64 {
    let offset0 = lidx.iter().rev().find(|&&v| v != 0).copied().unwrap_or(0);
    let bot = hts_bin_bot(bin) as usize;
    let val = if bot < lidx.len() { lidx[bot] } else { 0 };
    if val != 0 { val } else { offset0 }
}

// -----------------------------------------------------------------
// Binning scheme: htslib reg2bin – finest-level first
// -----------------------------------------------------------------

/// Compute the bin number for a 0-based half-open interval [beg, end).
fn reg2bin(beg: u64, end: u64) -> u32 {
    let e = end.saturating_sub(1);
    let mut s: u32 = MIN_SHIFT; // 14
    let mut t: u64 = ((1u64 << (3 * N_LVLS + 3)) - 1) / 7;
    for l in (1..=N_LVLS).rev() {
        t -= 1u64 << (3 * l);
        if (beg >> s) == (e >> s) {
            return (t + (beg >> s)) as u32;
        }
        s += 3;
    }
    0
}

// -----------------------------------------------------------------
// Index data structures
// -----------------------------------------------------------------

#[derive(Clone)]
struct Chunk {
    start: u64,
    end: u64,
}

struct SeqIdx {
    name: String,
    bins: HashMap<u32, Vec<Chunk>>,
    lidx: Vec<u64>,
    min_voff: u64,
    max_voff: u64,
    n_mapped: u64,
}

impl SeqIdx {
    fn new(name: String) -> Self {
        SeqIdx {
            name,
            bins: HashMap::new(),
            lidx: Vec::new(),
            min_voff: u64::MAX,
            max_voff: 0,
            n_mapped: 0,
        }
    }

    fn add_chunk(&mut self, bin: u32, chunk: Chunk) {
        if chunk.start < self.min_voff {
            self.min_voff = chunk.start;
        }
        if chunk.end > self.max_voff {
            self.max_voff = chunk.end;
        }
        self.n_mapped += 1;
        self.bins.entry(bin).or_default().push(chunk);
    }

    fn update_lidx(&mut self, beg: u64, end: u64, voff: u64) {
        if end == 0 {
            return;
        }
        let win_beg = (beg >> MIN_SHIFT) as usize;
        let win_end = ((end - 1) >> MIN_SHIFT) as usize;
        if win_end >= self.lidx.len() {
            self.lidx.resize(win_end + 1, 0);
        }
        for i in win_beg..=win_end {
            if self.lidx[i] == 0 {
                self.lidx[i] = voff;
            }
        }
    }
}

// -----------------------------------------------------------------
// compress_binning
// -----------------------------------------------------------------

fn merge_chunks_block_adjacent(chunks: &mut Vec<Chunk>) {
    if chunks.len() <= 1 {
        return;
    }
    chunks.sort_unstable_by_key(|c| c.start);
    let mut out: Vec<Chunk> = Vec::new();
    for c in chunks.drain(..) {
        match out.last_mut() {
            Some(last) if c.start <= last.end.saturating_add(HTS_MIN_MARKER_DIST) => {
                if c.end > last.end {
                    last.end = c.end;
                }
            }
            _ => out.push(c),
        }
    }
    *chunks = out;
}

fn compress_binning(bins: &mut HashMap<u32, Vec<Chunk>>) {
    for chunks in bins.values_mut() {
        chunks.sort_unstable_by_key(|c| c.start);
    }

    for l in (1..=N_LVLS).rev() {
        let level_first = hts_bin_first(l);
        let level_last = hts_bin_first(l + 1);

        let candidates: Vec<u32> = bins
            .keys()
            .filter(|&&b| b >= level_first && b < level_last)
            .cloned()
            .collect();

        for b in candidates {
            let parent = hts_bin_parent(b);
            if !bins.contains_key(&parent) {
                continue;
            }
            let chunks = bins.get(&b).unwrap();
            if chunks.is_empty() {
                continue;
            }
            let first_start = chunks.iter().map(|c| c.start).min().unwrap();
            let last_end   = chunks.iter().map(|c| c.end).max().unwrap();
            let span = (last_end >> 16).saturating_sub(first_start >> 16);
            if span < HTS_MIN_MARKER_DIST {
                let child_chunks = bins.remove(&b).unwrap();
                let parent_chunks = bins.get_mut(&parent).unwrap();
                parent_chunks.extend(child_chunks);
                parent_chunks.sort_unstable_by_key(|c| c.start);
            }
        }
    }

    for chunks in bins.values_mut() {
        merge_chunks_block_adjacent(chunks);
    }
}

// -----------------------------------------------------------------
// Public API
// -----------------------------------------------------------------

/// Build a CSI index for a BGZF-compressed GFF3 file.
///
/// Reads from `bgzf_input` (a BGZF-compressed byte stream) and writes the
/// binary `.csi` index to `csi_output`.
pub fn csi_index_gff<R: Read, W: Write>(bgzf_input: R, csi_output: W) -> io::Result<()> {
    let mut reader = BgzfReader::new(bgzf_input);

    let mut seqs: Vec<SeqIdx> = Vec::new();
    let mut seq_map: HashMap<String, usize> = HashMap::new();

    let mut line_buf = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let (n, voff_start) = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break;
        }

        let line = strip_newline(&line_buf);

        if line.is_empty() || line[0] == b'#' {
            continue;
        }

        let fields: Vec<&[u8]> = line.splitn(6, |&b| b == b'\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let seqname = std::str::from_utf8(fields[0])
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF8 sequence name"))?
            .to_owned();

        let start_1: u64 = parse_u64(fields[3])?;
        let end_1: u64 = parse_u64(fields[4])?;

        // GFF3 columns are 1-based, inclusive ⟹ convert to 0-based half-open
        let beg = start_1.saturating_sub(1);
        let end = end_1;

        let voff_end = reader.virtual_offset();
        let bin = reg2bin(beg, end);

        let tid = match seq_map.get(&seqname) {
            Some(&id) => id,
            None => {
                let id = seqs.len();
                seqs.push(SeqIdx::new(seqname.clone()));
                seq_map.insert(seqname, id);
                id
            }
        };

        let chunk = Chunk { start: voff_start, end: voff_end };
        seqs[tid].add_chunk(bin, chunk);
        seqs[tid].update_lidx(beg, end, voff_start);
    }

    // Fill trailing zeros in lidx.
    let eof_voff = reader.virtual_offset();
    for seq in &mut seqs {
        let mut seen_nonzero = false;
        for slot in seq.lidx.iter_mut() {
            if *slot != 0 {
                seen_nonzero = true;
            } else if seen_nonzero {
                *slot = eof_voff;
            }
        }
    }

    // Apply compress_binning and inject pseudo-bin META_BIN per sequence.
    for seq in &mut seqs {
        compress_binning(&mut seq.bins);

        let min_voff = if seq.min_voff == u64::MAX { 0 } else { seq.min_voff };
        seq.bins.insert(
            META_BIN,
            vec![
                Chunk { start: min_voff,    end: seq.max_voff },
                Chunk { start: seq.n_mapped, end: 0 },
            ],
        );
    }

    // -----------------------------------------------------------------
    // Write the .csi binary format (all little-endian), BGZF-compressed
    // -----------------------------------------------------------------
    let mut w = BgzfWriter::new(csi_output);

    // Magic
    w.write_all(b"CSI\x01")?;

    // min_shift, n_lvls
    w.write_all(&(MIN_SHIFT as i32).to_le_bytes())?;
    w.write_all(&(N_LVLS as i32).to_le_bytes())?;

    // Build names buffer for meta section
    let mut names_buf: Vec<u8> = Vec::new();
    for seq in &seqs {
        names_buf.extend_from_slice(seq.name.as_bytes());
        names_buf.push(0);
    }
    let l_nm = names_buf.len() as u32;

    // l_meta = 7 u32 fields (28 bytes) + names blob
    let l_meta: u32 = 28 + l_nm;
    w.write_all(&l_meta.to_le_bytes())?;

    // Meta blob: preset, col_seq, col_beg, col_end, meta_char, line_skip, l_nm, names.
    w.write_all(&0u32.to_le_bytes())?;   // preset = TBX_GENERIC
    w.write_all(&1u32.to_le_bytes())?;   // col_seq = 1 (1-based)
    w.write_all(&4u32.to_le_bytes())?;   // col_beg = 4 (1-based)
    w.write_all(&5u32.to_le_bytes())?;   // col_end = 5 (1-based)
    w.write_all(&35u32.to_le_bytes())?;  // meta_char = '#'
    w.write_all(&0u32.to_le_bytes())?;   // line_skip = 0
    w.write_all(&l_nm.to_le_bytes())?;   // l_nm
    w.write_all(&names_buf)?;            // seq names (null-terminated, concatenated)

    // n_ref
    w.write_all(&(seqs.len() as i32).to_le_bytes())?;

    // Per-sequence index data
    for seq in &seqs {
        let mut bin_ids: Vec<u32> = seq.bins.keys().cloned().collect();
        bin_ids.sort_unstable();

        w.write_all(&(bin_ids.len() as i32).to_le_bytes())?;
        for bin in &bin_ids {
            let chunks = &seq.bins[bin];
            let loff = compute_loff(*bin, &seq.lidx);
            w.write_all(&bin.to_le_bytes())?;
            w.write_all(&loff.to_le_bytes())?;
            w.write_all(&(chunks.len() as i32).to_le_bytes())?;
            for chunk in chunks {
                w.write_all(&chunk.start.to_le_bytes())?;
                w.write_all(&chunk.end.to_le_bytes())?;
            }
        }
    }

    // n_no_coor = 0
    w.write_all(&0u64.to_le_bytes())?;
    w.finish()?;

    Ok(())
}

// -----------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------

fn strip_newline(buf: &[u8]) -> &[u8] {
    let mut end = buf.len();
    while end > 0 && (buf[end - 1] == b'\n' || buf[end - 1] == b'\r') {
        end -= 1;
    }
    &buf[..end]
}

fn parse_u64(bytes: &[u8]) -> io::Result<u64> {
    let s = std::str::from_utf8(bytes)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "non-UTF8 field"))?
        .trim();
    s.parse::<u64>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("cannot parse integer: {}", s)))
}
