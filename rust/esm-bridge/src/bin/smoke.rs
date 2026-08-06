//! Native smoke test / benchmark: writes a `.npy` matrix, optionally against a reference, just embeddings

// Native-only: `wasm-pack` builds every target, so the whole body is cfg'd out.
#[cfg(not(target_family = "wasm"))]
mod native {

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use esm_bridge::{fasta, BatchPolicy, Embedder, HIDDEN};

#[derive(Parser, Debug)]
#[command(about = "Embed a protein FASTA with ESM-2 (8M) on the CPU backend")]
struct Args {
    /// Input protein FASTA.
    #[arg(long)]
    fasta: PathBuf,

    /// Burnpack weights (defaults to the build.rs byproduct).
    #[arg(long, default_value = "model/generated/esm2_t6_8M_UR50D.bpk")]
    weights: PathBuf,

    /// Where to write the [n, 320] float32 .npy matrix.
    #[arg(long)]
    out: Option<PathBuf>,

    /// Reference .npy to compare against.
    #[arg(long)]
    reference: Option<PathBuf>,

    /// Force a fixed batch size. Omit for the automatic batching the browser uses.
    #[arg(long)]
    batch_size: Option<usize>,

    #[arg(long, default_value_t = 1022)]
    max_len: usize,

    /// Re-run with batch size 1 and assert the result is identical (padding invariance).
    #[arg(long)]
    check_padding_invariance: bool,

    /// Which backend to run on. `webgpu` is the closest native proxy for the browser.
    #[arg(long, default_value = "cpu", value_parser = ["cpu", "webgpu"])]
    backend: String,
}

/// Minimal NumPy .npy v1.0 writer, little-endian float32, C order.
fn write_npy(path: &PathBuf, data: &[f32], shape: [usize; 2]) -> std::io::Result<()> {
    let dict = format!(
        "{{'descr': '<f4', 'fortran_order': False, 'shape': ({}, {}), }}",
        shape[0], shape[1]
    );
    // magic(6) + version(2) + headerLen(2) = 10; the total prefix must be a multiple of 64.
    let padding = (64 - ((10 + dict.len() + 1) % 64)) % 64;
    let header = format!("{dict}{}\n", " ".repeat(padding));

    let mut f = fs::File::create(path)?;
    f.write_all(b"\x93NUMPY")?;
    f.write_all(&[1u8, 0u8])?;
    f.write_all(&(header.len() as u16).to_le_bytes())?;
    f.write_all(header.as_bytes())?;
    for v in data {
        f.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

/// Minimal .npy reader for the exact subset we write.
fn read_npy(path: &PathBuf) -> Result<(Vec<f32>, [usize; 2]), String> {
    let buf = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if &buf[0..6] != b"\x93NUMPY" {
        return Err("not a .npy file".into());
    }
    let header_len = u16::from_le_bytes([buf[8], buf[9]]) as usize;
    let header = String::from_utf8_lossy(&buf[10..10 + header_len]).to_string();
    if !header.contains("'<f4'") {
        return Err(format!("expected float32 little-endian, got header: {header}"));
    }
    let shape_str = header
        .split("'shape':")
        .nth(1)
        .and_then(|s| s.split('(').nth(1))
        .and_then(|s| s.split(')').next())
        .ok_or("could not parse shape")?;
    let dims: Vec<usize> = shape_str
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    if dims.len() != 2 {
        return Err(format!("expected a 2-D array, got {dims:?}"));
    }
    let body = &buf[10 + header_len..];
    let vals: Vec<f32> = body
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();
    Ok((vals, [dims[0], dims[1]]))
}

fn cosine(a: &[f32], b: &[f32]) -> f64 {
    let (mut dot, mut na, mut nb) = (0f64, 0f64, 0f64);
    for (x, y) in a.iter().zip(b) {
        dot += (*x as f64) * (*y as f64);
        na += (*x as f64).powi(2);
        nb += (*y as f64).powi(2);
    }
    if na == 0.0 || nb == 0.0 {
        return 0.0;
    }
    dot / (na.sqrt() * nb.sqrt())
}

pub fn run() -> Result<(), String> {
    let args = Args::parse();

    let weights = fs::read(&args.weights)
        .map_err(|e| format!("{}: {e}\nRun `cargo build` first to generate it.", args.weights.display()))?;
    println!("weights   : {} ({:.1} MB)", args.weights.display(), weights.len() as f64 / 1048576.0);

    let fasta_bytes = fs::read(&args.fasta).map_err(|e| format!("{}: {e}", args.fasta.display()))?;
    let records = fasta::parse(&fasta_bytes)?;
    println!("sequences : {}", records.len());

    let t_load = Instant::now();
    let embedder = if args.backend == "webgpu" {
        let e = pollster::block_on(Embedder::new_gpu(weights, None, false));
        // Native Vulkan exposes SHADER_INT64 and browsers do not, so passing here does not
        // prove the browser will be happy — see build.rs.
        pollster::block_on(e.warmup_checked())?;
        e
    } else {
        Embedder::new_cpu(weights)
    };
    println!("backend   : {} (loaded in {:.2} s)", embedder.backend(), t_load.elapsed().as_secs_f64());

    let policy = match args.batch_size {
        Some(n) => BatchPolicy::Fixed(n),
        None => BatchPolicy::Auto,
    };
    println!(
        "budget    : {} M attention elements, batch ceiling {}",
        embedder.budget_elems() / 1_000_000,
        embedder.max_batch()
    );

    let t0 = Instant::now();
    let (vectors, meta) = pollster::block_on(embedder.embed_records(
        &records,
        args.max_len,
        policy,
        |done, total| {
            if done == total {
                println!("progress  : {done}/{total}");
            }
        },
    ));
    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "batching  : {:?}, groups of {}-{} sequences",
        policy, meta.batch_min, meta.batch_max
    );

    let residues: usize = meta.lengths.iter().sum();
    println!(
        "embedded  : {} seqs, {} residues in {:.2} s ({:.1} ms/seq, {:.0} residues/s)",
        meta.n,
        residues,
        elapsed,
        1000.0 * elapsed / meta.n as f64,
        residues as f64 / elapsed
    );
    if meta.truncated.iter().any(|&t| t) {
        let n = meta.truncated.iter().filter(|&&t| t).count();
        println!("truncated : {n} sequence(s) hit --max-len {}", args.max_len);
    }

    for (i, id) in meta.ids.iter().enumerate() {
        let row = &vectors[i * HIDDEN..(i + 1) * HIDDEN];
        let norm: f32 = row.iter().map(|v| v * v).sum::<f32>().sqrt();
        println!(
            "  {:<40} len={:<5} |v|={:.4} v[0..3]={:?}",
            id,
            meta.lengths[i],
            norm,
            &row[..3].iter().map(|v| (v * 1e4).round() / 1e4).collect::<Vec<_>>()
        );
        if !row.iter().all(|v| v.is_finite()) {
            return Err(format!("row {i} ({id}) contains non-finite values"));
        }
    }

    if let Some(out) = &args.out {
        write_npy(out, &vectors, [meta.n, HIDDEN]).map_err(|e| e.to_string())?;
        println!("wrote     : {} [{}, {}]", out.display(), meta.n, HIDDEN);
    }

    let mut failed = false;

    if args.check_padding_invariance {
        println!("\n--- padding invariance ({policy:?} vs batch 1) ---");
        let (single, _) = pollster::block_on(embedder.embed_records(
            &records,
            args.max_len,
            BatchPolicy::Fixed(1),
            |_, _| {},
        ));
        let mut worst = 0f32;
        let mut worst_id = String::new();
        for (i, id) in meta.ids.iter().enumerate() {
            let d = vectors[i * HIDDEN..(i + 1) * HIDDEN]
                .iter()
                .zip(&single[i * HIDDEN..(i + 1) * HIDDEN])
                .map(|(a, b)| (a - b).abs())
                .fold(0f32, f32::max);
            if d > worst {
                worst = d;
                worst_id = id.clone();
            }
        }
        println!("max abs diff : {worst:.3e} (worst: {worst_id})");
        if worst > 1e-5 {
            println!("FAIL: padding changes the embedding — attention/pool mask bug");
            failed = true;
        } else {
            println!("PASS");
        }
    }

    if let Some(refp) = &args.reference {
        println!("\n--- vs reference {} ---", refp.display());
        let (refv, shape) = read_npy(refp)?;
        if shape != [meta.n, HIDDEN] {
            return Err(format!("reference shape {shape:?} != [{}, {HIDDEN}]", meta.n));
        }
        let mut worst_abs = 0f64;
        let mut worst_cos = 1f64;
        for i in 0..meta.n {
            let a = &vectors[i * HIDDEN..(i + 1) * HIDDEN];
            let b = &refv[i * HIDDEN..(i + 1) * HIDDEN];
            let abs = a.iter().zip(b).map(|(x, y)| (x - y).abs() as f64).fold(0f64, f64::max);
            let cos = cosine(a, b);
            println!("  {:<40} maxabs={:.3e}  cos={:.8}", meta.ids[i], abs, cos);
            worst_abs = worst_abs.max(abs);
            worst_cos = worst_cos.min(cos);
        }
        println!("worst max-abs-diff : {worst_abs:.3e}  (threshold 1e-3)");
        println!("worst cosine       : {worst_cos:.8}  (threshold 0.9999)");
        if worst_abs >= 1e-3 || worst_cos <= 0.9999 {
            println!("FAIL");
            failed = true;
        } else {
            println!("PASS");
        }
    }

    if failed {
        return Err("one or more checks failed".into());
    }
    Ok(())
}

} // mod native

#[cfg(not(target_family = "wasm"))]
fn main() -> Result<(), String> {
    native::run()
}

#[cfg(target_family = "wasm")]
fn main() {}
