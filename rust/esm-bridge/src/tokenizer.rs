//! Hand-written ESM-2 tokenizer.

// i32, not i64: WGSL has no 64-bit integer type. See the note in build.rs.
pub const CLS: i32 = 0;
pub const PAD: i32 = 1;
pub const EOS: i32 = 2;
pub const UNK: i32 = 3;
pub const MASK: i32 = 32;
pub const VOCAB_SIZE: usize = 33;

/// Index == token id. Order is prepend_toks + standard_toks + append_toks.
pub const VOCAB: [&str; VOCAB_SIZE] = [
    "<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E", "R", "T", "I", "D", "P",
    "K", "Q", "N", "F", "Y", "M", "H", "W", "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>",
    "<mask>",
];

/// `max_position_embeddings` is 1026, which leaves 1024 tokens, i.e. 1022 residues once
/// `<cls>` and `<eos>` are accounted for.
pub const MAX_RESIDUES: usize = 1022;

/// ASCII byte -> token id; 255 marks "not in the alphabet" (mapped to [`UNK`]).
const fn build_lut() -> [u8; 256] {
    let mut l = [255u8; 256];
    l[b'L' as usize] = 4;
    l[b'A' as usize] = 5;
    l[b'G' as usize] = 6;
    l[b'V' as usize] = 7;
    l[b'S' as usize] = 8;
    l[b'E' as usize] = 9;
    l[b'R' as usize] = 10;
    l[b'T' as usize] = 11;
    l[b'I' as usize] = 12;
    l[b'D' as usize] = 13;
    l[b'P' as usize] = 14;
    l[b'K' as usize] = 15;
    l[b'Q' as usize] = 16;
    l[b'N' as usize] = 17;
    l[b'F' as usize] = 18;
    l[b'Y' as usize] = 19;
    l[b'M' as usize] = 20;
    l[b'H' as usize] = 21;
    l[b'W' as usize] = 22;
    l[b'C' as usize] = 23;
    l[b'X' as usize] = 24;
    l[b'B' as usize] = 25;
    l[b'U' as usize] = 26;
    l[b'Z' as usize] = 27;
    l[b'O' as usize] = 28;
    l[b'.' as usize] = 29;
    l[b'-' as usize] = 30;
    l
}

pub static AA_LUT: [u8; 256] = build_lut();

/// Iterator over the residues of a raw FASTA sequence: whitespace is skipped and a single
/// trailing `*` (stop codon, as emitted by prodigal/orphos) is dropped.
fn residues(seq: &[u8]) -> impl Iterator<Item = u8> + '_ {
    let last_non_ws = seq.iter().rposition(|c| !c.is_ascii_whitespace());
    let drop_star = matches!(last_non_ws, Some(i) if seq[i] == b'*');
    let end = if drop_star { last_non_ws.unwrap() } else { seq.len() };
    seq[..end]
        .iter()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| c.to_ascii_uppercase())
}

/// Number of residues that will actually be tokenised.
pub fn residue_count(seq: &[u8]) -> usize {
    residues(seq).count()
}

#[derive(Debug, Clone)]
pub struct Encoded {
    /// Row-major `[batch, len]`.
    pub input_ids: Vec<i32>,
    /// 1 on `<cls>`, residues and `<eos>`; 0 on `<pad>`.
    pub attention_mask: Vec<i32>,
    /// 1.0 on residues only (excludes `<cls>`/`<eos>`/`<pad>`) — used for mean pooling.
    pub pool_mask: Vec<f32>,
    pub batch: usize,
    pub len: usize,
    /// Whether each sequence was longer than `max_residues` and got head-truncated.
    pub truncated: Vec<bool>,
    /// Residues actually kept per sequence.
    pub kept: Vec<usize>,
}

/// Encode a batch, padding to the longest member.
pub fn encode_batch(seqs: &[&[u8]], max_residues: usize) -> Encoded {
    let batch = seqs.len();
    let kept: Vec<usize> = seqs
        .iter()
        .map(|s| residue_count(s).min(max_residues))
        .collect(); // how many we reduce
    let max_kept = kept.iter().copied().max().unwrap_or(0);
    let len = max_kept + 2; // <cls> ... <eos>

    let mut input_ids = vec![PAD; batch * len];
    let mut attention_mask = vec![0i32; batch * len];
    let mut pool_mask = vec![0f32; batch * len];
    let mut truncated = vec![false; batch];

    for (b, seq) in seqs.iter().enumerate() {
        let off = b * len;
        input_ids[off] = CLS;
        attention_mask[off] = 1;

        let mut k = 1usize;
        for c in residues(seq).take(max_residues) {
            let id = AA_LUT[c as usize];
            input_ids[off + k] = if id == 255 { UNK } else { id as i32 };
            attention_mask[off + k] = 1;
            pool_mask[off + k] = 1.0;
            k += 1;
        }
        input_ids[off + k] = EOS;
        attention_mask[off + k] = 1;

        truncated[b] = residue_count(seq) > max_residues;
    }

    Encoded {
        input_ids,
        attention_mask,
        pool_mask,
        batch,
        len,
        truncated,
        kept,
    }
}

/// Round a padded length up to a multiple of this.
pub const LEN_QUANTUM: usize = 128;

/// Padded length actually used for a batch whose longest member has `max_residues` residues.
pub fn padded_len(max_residues: usize) -> usize {
    let needed = max_residues + 2; // <cls> .. <eos>
    needed.div_ceil(LEN_QUANTUM) * LEN_QUANTUM
}

/// Largest batch whose attention scores fit the budget at a given padded length.
pub fn batch_cap(padded: usize, budget_elems: usize, heads: usize, max_batch: usize) -> usize {
    (budget_elems / (heads * padded * padded)).clamp(1, max_batch)
}

/// Length-sorted bucketing with an automatically chosen batch size.
pub fn plan_batches_auto(
    residue_counts: &[usize],
    budget_elems: usize,
    heads: usize,
    max_batch: usize,
) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..residue_counts.len()).collect();
    order.sort_unstable_by_key(|&i| residue_counts[i]);

    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut cur: Vec<usize> = Vec::new();
    let mut cur_bin = 0usize;

    for &i in &order {
        let bin = padded_len(residue_counts[i]);
        let cap = batch_cap(bin, budget_elems, heads, max_batch);
        if !cur.is_empty() && (bin != cur_bin || cur.len() + 1 > cap) {
            out.push(std::mem::take(&mut cur));
        }
        cur_bin = bin;
        cur.push(i);
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Fixed-size length-sorted bucketing, for tests and benchmarks only.
pub fn plan_batches_fixed(residue_counts: &[usize], batch_size: usize) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..residue_counts.len()).collect();
    order.sort_unstable_by_key(|&i| residue_counts[i]);
    order
        .chunks(batch_size.max(1))
        .map(|c| c.to_vec())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids_of(seq: &[u8]) -> Vec<i32> {
        let e = encode_batch(&[seq], MAX_RESIDUES);
        e.input_ids
    }

    #[test]
    fn vocab_is_well_formed() {
        assert_eq!(VOCAB.len(), VOCAB_SIZE);
        assert_eq!(VOCAB[CLS as usize], "<cls>");
        assert_eq!(VOCAB[PAD as usize], "<pad>");
        assert_eq!(VOCAB[EOS as usize], "<eos>");
        assert_eq!(VOCAB[UNK as usize], "<unk>");
        assert_eq!(VOCAB[MASK as usize], "<mask>");
        // Every single-character vocab entry must round-trip through the LUT.
        for (i, tok) in VOCAB.iter().enumerate() {
            if tok.len() == 1 {
                let c = tok.as_bytes()[0];
                assert_eq!(AA_LUT[c as usize] as usize, i, "LUT mismatch for {tok}");
            }
        }
    }

    #[test]
    fn wraps_with_cls_and_eos() {
        // M=20, K=15
        assert_eq!(ids_of(b"MK"), vec![CLS, 20, 15, EOS]);
    }

    #[test]
    fn lowercase_is_uppercased_not_unk() {
        assert_eq!(ids_of(b"mk"), ids_of(b"MK"));
    }

    #[test]
    fn rare_residues_keep_their_own_ids() {
        // X=24, B=25, U=26, Z=27, O=28 — must NOT collapse to X or <unk>.
        assert_eq!(ids_of(b"XBUZO"), vec![CLS, 24, 25, 26, 27, 28, EOS]);
    }

    #[test]
    fn unknown_characters_become_unk() {
        // 'J' and '1' are not in the ESM alphabet.
        assert_eq!(ids_of(b"J1"), vec![CLS, UNK, UNK, EOS]);
    }

    #[test]
    fn trailing_stop_codon_is_dropped_but_internal_star_is_unk() {
        assert_eq!(ids_of(b"MK*"), vec![CLS, 20, 15, EOS]);
        assert_eq!(ids_of(b"MK*\n"), vec![CLS, 20, 15, EOS]);
        // Only ONE trailing star is dropped.
        assert_eq!(ids_of(b"MK**"), vec![CLS, 20, 15, UNK, EOS]);
        // A star in the middle is a genuine unknown.
        assert_eq!(ids_of(b"M*K"), vec![CLS, 20, UNK, 15, EOS]);
    }

    #[test]
    fn whitespace_and_wrapping_are_ignored() {
        assert_eq!(ids_of(b"M K\nM\r\nK"), ids_of(b"MKMK"));
    }

    #[test]
    fn padding_and_masks_are_consistent() {
        let e = encode_batch(&[b"MK".as_slice(), b"MKMKM".as_slice()], MAX_RESIDUES);
        assert_eq!(e.batch, 2);
        assert_eq!(e.len, 7); // 5 residues + cls + eos
        assert_eq!(e.kept, vec![2, 5]);

        // Row 0: <cls> M K <eos> <pad> <pad> <pad>
        assert_eq!(&e.input_ids[0..7], &[CLS, 20, 15, EOS, PAD, PAD, PAD]);
        assert_eq!(&e.attention_mask[0..7], &[1, 1, 1, 1, 0, 0, 0]);
        assert_eq!(&e.pool_mask[0..7], &[0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]);

        // pool_mask must count exactly the residues, never cls/eos/pad.
        for b in 0..e.batch {
            let row = &e.pool_mask[b * e.len..(b + 1) * e.len];
            assert_eq!(row.iter().sum::<f32>() as usize, e.kept[b]);
        }
        // attention_mask is always pool_mask + 2 (cls and eos).
        for b in 0..e.batch {
            let row = &e.attention_mask[b * e.len..(b + 1) * e.len];
            assert_eq!(row.iter().sum::<i32>() as usize, e.kept[b] + 2);
        }
    }

    #[test]
    fn truncation_is_head_truncation_and_is_flagged() {
        let seq = vec![b'A'; 10];
        let e = encode_batch(&[seq.as_slice()], 4);
        assert_eq!(e.len, 6);
        assert_eq!(e.kept, vec![4]);
        assert_eq!(e.truncated, vec![true]);
        assert_eq!(e.input_ids, vec![CLS, 5, 5, 5, 5, EOS]);

        let e2 = encode_batch(&[b"AAAA".as_slice()], 4);
        assert_eq!(e2.truncated, vec![false]);
    }

    #[test]
    fn empty_sequence_is_just_cls_eos() {
        let e = encode_batch(&[b"".as_slice()], MAX_RESIDUES);
        assert_eq!(e.input_ids, vec![CLS, EOS]);
        assert_eq!(e.kept, vec![0]);
    }

    #[test]
    fn plan_batches_fixed_sorts_by_length_and_covers_every_index() {
        let counts = [100, 5, 50, 5, 900];
        let plan = plan_batches_fixed(&counts, 2);
        let mut seen: Vec<usize> = plan.iter().flatten().copied().collect();
        seen.sort_unstable();
        assert_eq!(seen, vec![0, 1, 2, 3, 4]);
        // Shortest sequences must land in the first bucket.
        let first: Vec<usize> = plan[0].iter().map(|&i| counts[i]).collect();
        assert_eq!(first, vec![5, 5]);
    }

    // A fixed budget, so these pin batch_cap's arithmetic rather than the production
    // constants, which now differ per backend and per target.
    const BUDGET: usize = 16_000_000;
    const HEADS: usize = 20;

    #[test]
    fn padded_len_quantises_upwards() {
        assert_eq!(padded_len(1), LEN_QUANTUM); // 3 tokens -> 128
        assert_eq!(padded_len(126), 128); // exactly 128 tokens
        assert_eq!(padded_len(127), 256); // 129 tokens -> next quantum
        assert_eq!(padded_len(1022), 1024); // the model maximum lands exactly
        // Never truncates: the padded length always covers the tokens needed.
        for r in [0, 1, 63, 200, 511, 1022] {
            assert!(padded_len(r) >= r + 2, "padded_len({r}) too small");
        }
    }

    #[test]
    fn batch_cap_falls_with_the_square_of_length() {
        // 20 heads * 1024^2 * 1 = 20.9M elements > 16M, so the maximum length forces 1.
        assert_eq!(batch_cap(1024, BUDGET, HEADS, 32), 1);
        assert_eq!(batch_cap(512, BUDGET, HEADS, 32), 3);
        assert_eq!(batch_cap(256, BUDGET, HEADS, 32), 12);
        // Short sequences are limited by max_batch, not by memory.
        assert_eq!(batch_cap(128, BUDGET, HEADS, 32), 32);
        // ...and the ceiling is honoured whatever the budget.
        assert_eq!(batch_cap(128, usize::MAX / 64, HEADS, 8), 8);
    }

    #[test]
    fn auto_plan_covers_every_index_exactly_once() {
        let counts: Vec<usize> = (0..101).map(|i| 38 + (i * 13) % 1400).collect();
        let plan = plan_batches_auto(&counts, BUDGET, HEADS, 32);
        let mut seen: Vec<usize> = plan.iter().flatten().copied().collect();
        seen.sort_unstable();
        assert_eq!(seen, (0..101).collect::<Vec<_>>());
        assert!(plan.iter().all(|g| !g.is_empty()), "no empty groups");
    }

    #[test]
    fn auto_plan_never_exceeds_the_budget_or_the_ceiling() {
        let counts: Vec<usize> = (0..200).map(|i| 10 + (i * 37) % 1013).collect();
        let max_batch = 32;
        for group in plan_batches_auto(&counts, BUDGET, HEADS, max_batch) {
            let padded = padded_len(group.iter().map(|&i| counts[i]).max().unwrap());
            let elems = group.len() * HEADS * padded * padded;
            assert!(
                group.len() <= max_batch,
                "group of {} exceeds the ceiling",
                group.len()
            );
            // A single sequence is always allowed through even if it alone exceeds the
            // budget, since there is no smaller unit of work.
            assert!(
                group.len() == 1 || elems <= BUDGET,
                "group of {} at padded {padded} needs {elems} elements, over budget",
                group.len()
            );
        }
    }

    #[test]
    fn auto_plan_batches_long_sequences_alone_and_short_ones_together() {
        let long = vec![1022usize; 4];
        assert!(plan_batches_auto(&long, BUDGET, HEADS, 32)
            .iter()
            .all(|g| g.len() == 1));

        let short = vec![60usize; 100];
        let plan = plan_batches_auto(&short, BUDGET, HEADS, 32);
        // 60 residues -> padded 128 -> cap is the ceiling, so groups fill to 32.
        assert_eq!(plan[0].len(), 32);
        assert_eq!(plan.len(), 4); // 32 + 32 + 32 + 4
    }

    #[test]
    fn auto_plan_groups_share_one_padded_length() {
        // Wide spread on purpose: without binning, short sequences would land in a group
        // padded to a far longer member and waste most of their compute.
        let counts: Vec<usize> = (0..300).map(|i| 20 + (i * 7) % 1000).collect();
        for group in plan_batches_auto(&counts, BUDGET, HEADS, 32) {
            let lens: Vec<usize> = group.iter().map(|&i| padded_len(counts[i])).collect();
            assert!(
                lens.windows(2).all(|w| w[0] == w[1]),
                "group mixes padded lengths {lens:?} — padding waste is unbounded"
            );
        }
    }

    #[test]
    fn auto_plan_handles_edge_cases() {
        assert!(plan_batches_auto(&[], BUDGET, HEADS, 32).is_empty());
        assert_eq!(plan_batches_auto(&[0], BUDGET, HEADS, 32), vec![vec![0]]);
        // A ceiling of 1 must give one group per sequence.
        assert_eq!(plan_batches_auto(&[10, 20, 30], BUDGET, HEADS, 1).len(), 3);
    }
}
