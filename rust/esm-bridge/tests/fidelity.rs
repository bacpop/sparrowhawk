//! Numerical fidelity/correctness test

use std::path::Path;

use esm_bridge::{fasta, BatchPolicy, Embedder, HIDDEN};

/// Serialises the GPU tests.
static GPU_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

const WEIGHTS: &str = "model/generated/esm2_t6_8M_UR50D.bpk";
const FASTA: &str = "tests/data/seqs.faa";
const REFERENCE: &str = "tests/data/reference.npy";

/// 20 proteins, 38–1463 residues, one over the 1022 cap. Small because `cargo test` builds
/// NdArray unoptimised.
const FASTA_20: &str = "tests/data/proteins20.faa";
/// Produced by `scripts/reference_official.py` from facebook/esm2_t6_8M_UR50D itself.
const REFERENCE_OFFICIAL: &str = "tests/data/reference_official.npy";


/// This magic reads Numpy .npy files
fn read_npy_f32(path: &str) -> (Vec<f32>, [usize; 2]) {
    let buf = std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    assert_eq!(&buf[0..6], b"\x93NUMPY", "{path} is not a .npy file");
    let header_len = u16::from_le_bytes([buf[8], buf[9]]) as usize;
    let header = String::from_utf8_lossy(&buf[10..10 + header_len]).to_string();
    assert!(header.contains("'<f4'"), "expected float32 LE: {header}");
    let shape_str = header
        .split("'shape':")
        .nth(1)
        .and_then(|s| s.split('(').nth(1))
        .and_then(|s| s.split(')').next())
        .expect("shape in npy header");
    let dims: Vec<usize> = shape_str
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();
    let vals: Vec<f32> = buf[10 + header_len..]
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    (vals, [dims[0], dims[1]])
}


fn embedder() -> Option<Embedder> {
    if !Path::new(WEIGHTS).exists() {
        eprintln!("skipping: {WEIGHTS} missing (run `cargo build` first)");
        return None;
    }
    Some(Embedder::new_cpu(std::fs::read(WEIGHTS).unwrap()))
}

/// Worst per-row absolute difference and worst per-row cosine, against a reference matrix.
/// Panics if any of our rows is non-finite, which is a different failure from being wrong.
fn worst_abs_and_cos(ours: &[f32], reference: &[f32], ids: &[String]) -> (f64, f64) {
    let (mut worst_abs, mut worst_cos) = (0f64, 1f64);
    for (i, id) in ids.iter().enumerate() {
        let a = &ours[i * HIDDEN..(i + 1) * HIDDEN];
        let b = &reference[i * HIDDEN..(i + 1) * HIDDEN];
        assert!(
            a.iter().all(|v| v.is_finite()),
            "row {i} ({id}) has non-finite values"
        );
        worst_abs = worst_abs.max(
            a.iter()
                .zip(b)
                .map(|(x, y)| (x - y).abs() as f64)
                .fold(0f64, f64::max),
        );
        let (mut dot, mut na, mut nb) = (0f64, 0f64, 0f64);
        for (x, y) in a.iter().zip(b) {
            dot += (*x as f64) * (*y as f64);
            na += (*x as f64).powi(2);
            nb += (*y as f64).powi(2);
        }
        worst_cos = worst_cos.min(dot / (na.sqrt() * nb.sqrt()));
    }
    (worst_abs, worst_cos)
}

#[test]
fn matches_onnxruntime_reference() {
    let Some(e) = embedder() else { return };
    let records = fasta::parse(&std::fs::read(FASTA).unwrap()).unwrap();
    let (ours, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));

    let (reference, shape) = read_npy_f32(REFERENCE);
    assert_eq!(shape, [meta.n, HIDDEN], "reference shape mismatch");
    let (worst_abs, worst_cos) = worst_abs_and_cos(&ours, &reference, &meta.ids);

    assert!(
        worst_abs < 1e-3,
        "max abs diff {worst_abs:.3e} >= 1e-3 — the ONNX graph is being miscompiled. \
         Check that build.rs still sets .simplify(false) (burn-onnx#436)."
    );
    assert!(
        worst_cos > 0.9999,
        "worst cosine {worst_cos:.8} <= 0.9999 — see burn-onnx#436"
    );
}

/// Against Meta's *weights*, not our own graph, which shares them with the model under test.
/// Regenerate the fixture with `scripts/reference_official.py` when the weights change.
#[test]
fn matches_official_meta_weights() {
    let Some(e) = embedder() else { return };
    let records = fasta::parse(&std::fs::read(FASTA_20).unwrap()).unwrap();
    // Auto, not Fixed: fixed groups pad to their longest member, which the clamp then splits
    // back into batches of one. `padding_does_not_change_embeddings` pins batch-invariance.
    let (ours, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Auto, |_, _| {}));

    let (reference, shape) = read_npy_f32(REFERENCE_OFFICIAL);
    assert_eq!(shape, [meta.n, HIDDEN], "reference shape mismatch");
    let (worst_abs, worst_cos) = worst_abs_and_cos(&ours, &reference, &meta.ids);
    println!("vs facebook/esm2_t6_8M_UR50D over {} proteins: max abs {worst_abs:.3e}, worst cosine {worst_cos:.8}", meta.n);

    // Thresholds set from the measured values, with headroom for f32 accumulation order.
    assert!(
        worst_abs < 1e-3,
        "max abs diff {worst_abs:.3e} >= 1e-3 against Meta's own weights — the export in \
         model/ does not faithfully reproduce facebook/esm2_t6_8M_UR50D."
    );
    assert!(
        worst_cos > 0.9999,
        "worst cosine {worst_cos:.8} <= 0.9999 against Meta's own weights"
    );
}

#[test]
fn padding_does_not_change_embeddings() {
    let Some(e) = embedder() else { return };
    // Deliberately length-mixed, so batching forces heavy padding on the short sequences.
    let records = fasta::parse(&std::fs::read(FASTA).unwrap()).unwrap();

    let (batched, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));
    let (single, _) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(1), |_, _| {}));

    for i in 0..meta.n {
        let d = batched[i * HIDDEN..(i + 1) * HIDDEN]
            .iter()
            .zip(&single[i * HIDDEN..(i + 1) * HIDDEN])
            .map(|(a, b)| (a - b).abs())
            .fold(0f32, f32::max);
        assert!(
            d < 1e-5,
            "sequence {} ({}) changed by {d:.3e} when padded — attention/pool mask bug",
            i,
            meta.ids[i]
        );
    }
}

#[test]
fn row_order_follows_input_order() {
    let Some(e) = embedder() else { return };
    let records = fasta::parse(&std::fs::read(FASTA).unwrap()).unwrap();
    let (straight, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));

    // Reverse the file: each row must follow its own record, despite length-sorted batching.
    let reversed: Vec<fasta::Record> = records.iter().rev().cloned().collect();
    let (rev_out, rev_meta) = pollster::block_on(e.embed_records(&reversed, 1022, BatchPolicy::Fixed(8), |_, _| {}));

    for (i, id) in meta.ids.iter().enumerate() {
        let j = rev_meta.ids.iter().position(|x| x == id).unwrap();
        let d = straight[i * HIDDEN..(i + 1) * HIDDEN]
            .iter()
            .zip(&rev_out[j * HIDDEN..(j + 1) * HIDDEN])
            .map(|(a, b)| (a - b).abs())
            .fold(0f32, f32::max);
        assert!(d < 1e-5, "row for {id} moved or changed ({d:.3e})");
    }
}

/// GPU counterpart of [`matches_onnxruntime_reference`].
#[test]
#[ignore = "requires a GPU"]
fn gpu_matches_onnxruntime_reference() {
    let _guard = GPU_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    if !Path::new(WEIGHTS).exists() {
        eprintln!("skipping: {WEIGHTS} missing (run `cargo build` first)");
        return;
    }
    let e = pollster::block_on(Embedder::new_gpu(std::fs::read(WEIGHTS).unwrap(), None, false));
    assert_eq!(e.backend(), "webgpu");
    // The self-check must not produce false positives on a backend that works, otherwise it
    // would silently demote every GPU user to the CPU.
    pollster::block_on(e.warmup_checked())
        .expect("self-check must pass on a working GPU");

    let records = fasta::parse(&std::fs::read(FASTA).unwrap()).unwrap();
    let (ours, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));

    let (reference, shape) = read_npy_f32(REFERENCE);
    assert_eq!(shape, [meta.n, HIDDEN]);

    let mut worst_abs = 0f64;
    let mut worst_cos = 1f64;
    for i in 0..meta.n {
        let a = &ours[i * HIDDEN..(i + 1) * HIDDEN];
        let b = &reference[i * HIDDEN..(i + 1) * HIDDEN];
        assert!(a.iter().all(|v| v.is_finite()), "row {i} has non-finite values");
        worst_abs = worst_abs.max(
            a.iter().zip(b).map(|(x, y)| (x - y).abs() as f64).fold(0f64, f64::max),
        );
        let (mut dot, mut na, mut nb) = (0f64, 0f64, 0f64);
        for (x, y) in a.iter().zip(b) {
            dot += (*x as f64) * (*y as f64);
            na += (*x as f64).powi(2);
            nb += (*y as f64).powi(2);
        }
        worst_cos = worst_cos.min(dot / (na.sqrt() * nb.sqrt()));
    }
    eprintln!("gpu vs reference: max abs diff {worst_abs:.3e}, worst cosine {worst_cos:.8}");
    assert!(worst_abs < 1e-3, "gpu max abs diff {worst_abs:.3e} >= 1e-3");
    assert!(worst_cos > 0.9999, "gpu worst cosine {worst_cos:.8} <= 0.9999");
}

/// The GPU and CPU backends must agree with each other, not merely each with the reference.
#[test]
#[ignore = "requires a GPU"]
fn gpu_and_cpu_agree() {
    let _guard = GPU_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    if !Path::new(WEIGHTS).exists() {
        return;
    }
    let weights = std::fs::read(WEIGHTS).unwrap();
    let records = fasta::parse(&std::fs::read(FASTA).unwrap()).unwrap();

    let cpu = Embedder::new_cpu(weights.clone());
    let gpu = pollster::block_on(Embedder::new_gpu(weights, None, false));

    let (c, meta) = pollster::block_on(cpu.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));
    let (g, _) = pollster::block_on(gpu.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));

    for i in 0..meta.n {
        let d = c[i * HIDDEN..(i + 1) * HIDDEN]
            .iter()
            .zip(&g[i * HIDDEN..(i + 1) * HIDDEN])
            .map(|(a, b)| (a - b).abs())
            .fold(0f32, f32::max);
        assert!(d < 1e-3, "backends disagree on {} by {d:.3e}", meta.ids[i]);
    }
}

/// One max-length attention matrix is 20 * 1024^2 * 4 B = ~84 MB, so two would exceed
/// WebGPU's 128 MiB maxStorageBufferBindingSize. This asserts
#[test]
#[ignore = "requires a GPU"]
fn gpu_handles_max_length_sequences() {
    let _guard = GPU_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    if !Path::new(WEIGHTS).exists() {
        return;
    }
    let e = pollster::block_on(Embedder::new_gpu(std::fs::read(WEIGHTS).unwrap(), None, false));
    let records: Vec<fasta::Record> = (0..4)
        .map(|i| fasta::Record {
            id: format!("long_{i}"),
            desc: String::new(),
            seq: vec![b'A'; 1022],
        })
        .collect();
    // Ask for batch 8 deliberately: the clamp must reduce it rather than fail.
    let (v, meta) = pollster::block_on(e.embed_records(&records, 1022, BatchPolicy::Fixed(8), |_, _| {}));
    assert_eq!(meta.lengths, vec![1022; 4]);
    assert!(meta.truncated.iter().all(|&t| !t));
    assert!(v.iter().all(|x| x.is_finite()), "non-finite values at max length");
}

/// The self-check must pass on a backend that is known good, otherwise it would fire
/// spuriously and push everyone onto the CPU.
#[test]
fn warmup_self_check_passes_on_cpu() {
    let Some(e) = embedder() else { return };
    pollster::block_on(e.warmup_checked())
        .expect("CPU self-check must pass; WARMUP_EXPECTED has drifted from the model");
}

/// Guards the baked-in constant against a silent drift if the model or pooling changes.
#[test]
fn warmup_constant_matches_the_model() {
    let Some(e) = embedder() else { return };
    let rec = fasta::Record {
        id: "warmup".into(),
        desc: String::new(),
        seq: b"MK".to_vec(),
    };
    let (v, _) = pollster::block_on(e.embed_records(&[rec], 1022, BatchPolicy::Fixed(1), |_, _| {}));
    // Same first-8-dims the self-check compares, at a far tighter tolerance than the
    // check itself uses.
    let expected = [
        0.076_747_07f32, -0.184_892_62, 0.428_570_5, 0.382_838_2,
        -0.132_557_48, -0.145_408_18, -0.344_398_62, -0.148_791_76,
    ];
    for (i, (&got, &want)) in v.iter().zip(expected.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-5,
            "dim {i}: {got:.8} vs {want:.8} — update WARMUP_EXPECTED in src/lib.rs"
        );
    }
}

#[test]
fn truncation_is_flagged_and_bounded() {
    let Some(e) = embedder() else { return };
    let long = fasta::Record {
        id: "long".into(),
        desc: String::new(),
        seq: vec![b'A'; 40],
    };
    let (_, meta) = pollster::block_on(e.embed_records(&[long], 16, BatchPolicy::Fixed(1), |_, _| {}));
    assert_eq!(meta.lengths, vec![16]);
    assert_eq!(meta.truncated, vec![true]);
}
