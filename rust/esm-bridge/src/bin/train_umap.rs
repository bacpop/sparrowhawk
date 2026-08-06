//! Build the 2-D encoder the web tab uses. That is a simple neural network of two hidden layers,
//! a multilayer perceptron. First, we need to embed the dataset to be used for the training, then,
//! fit the projection/encoder/NN/MLP, and finally we can check for potential overfitting.

// Native-only, like smoke.rs: wasm-pack builds every target in the crate.
#[cfg(not(target_family = "wasm"))]
mod native {

use std::collections::HashSet;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use burn::backend::Autodiff;
use burn::module::AutodiffModule;
use burn::nn::loss::{MseLoss, Reduction};
use burn::optim::{AdamConfig, GradientsParams, Optimizer};
use burn::prelude::*;
use burn_store::{BurnpackStore, ModuleSnapshot};
use clap::{Parser, Subcommand};
use esm_bridge::projector::{UmapMlp, N_DIM};
use esm_bridge::{fasta, BatchPolicy, Embedder, HIDDEN};

#[derive(Parser, Debug)]
#[command(about = "Build the 2-D UMAP encoder from a corpus of protein FASTAs")]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Embed a corpus of FASTAs into `<work-dir>/embeddings.f32`.
    Embed(EmbedArgs),
    /// Fit the encoder to the coordinates from `scripts/umap_coords.py`.
    Fit(FitArgs),
    /// Score an existing encoder on the same split `fit` would use, without retraining.
    Eval(EvalArgs),
}

#[derive(Parser, Debug)]
struct EvalArgs {
    #[arg(long, default_value = "model/generated/umap_corpus")]
    work_dir: PathBuf,

    #[arg(long, default_value = "model/generated/esm2_umap_encoder.bpk")]
    encoder: PathBuf,

    /// Must match the run that produced the encoder, or rows it trained on land in the test
    /// set and the score flatters it.
    #[arg(long, default_value_t = 0.15)] // by default, 15% test, 85% training
    holdout: f32,

    #[arg(long, default_value_t = 42)] // for reproducibility
    seed: u64,
}

#[derive(Parser, Debug)]
struct EmbedArgs {
    /// TSV of FASTA paths: `<path>[\t<label>]`, one per line, `#` for comments.
    /// Paths resolve relative to the TSV's own directory.
    #[arg(long)]
    fasta_list: PathBuf,

    /// Where the embedding matrix and its manifest live. Reused by --resume.
    #[arg(long, default_value = "model/generated/umap_corpus")]
    work_dir: PathBuf,

    #[arg(long, default_value = "model/generated/esm2_t6_8M_UR50D.bpk")]
    weights: PathBuf,

    /// Same truncation the browser applies, so training and inference see identical input.
    #[arg(long, default_value_t = 1022)]
    max_len: usize,

    /// Skip FASTAs already recorded in the manifest.
    #[arg(long)]
    resume: bool,

    /// Keep exact-duplicate sequences. They inflate local density and distort the layout.
    #[arg(long)]
    keep_duplicates: bool,

    /// Embed on the CPU rather than the GPU.
    #[arg(long)]
    cpu: bool,
}

#[derive(Parser, Debug)]
struct FitArgs {
    #[arg(long, default_value = "model/generated/umap_corpus")]
    work_dir: PathBuf,

    #[arg(long, default_value = "model/generated/esm2_umap_encoder.bpk")]
    out: PathBuf,

    #[arg(long, default_value_t = 500)]
    epochs: usize,

    #[arg(long, default_value_t = 256)]
    batch_size: usize,

    #[arg(long, default_value_t = 1e-3)]
    lr: f64,

    /// Test fraction (one minus gives you the train one)
    #[arg(long, default_value_t = 0.15)]
    holdout: f32,

    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Fit on the CPU rather than the GPU.
    #[arg(long)]
    cpu: bool,
}


// struct for reading the input data set (as fasta file with aminoacid sequences)
struct CorpusEntry {
    path: PathBuf,
    label: String,
}

fn read_corpus(tsv: &Path) -> Result<Vec<CorpusEntry>, String> {
    let base = tsv.parent().unwrap_or_else(|| Path::new("."));
    let text = fs::read_to_string(tsv).map_err(|e| format!("{}: {e}", tsv.display()))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.split('\t');
        let path = cols.next().unwrap_or_default().trim();
        if path.is_empty() {
            continue;
        }
        out.push(CorpusEntry {
            path: base.join(path),
            label: cols.next().unwrap_or_default().trim().to_string(),
        });
    }
    if out.is_empty() {
        return Err(format!("{} lists no FASTA files", tsv.display()));
    }
    Ok(out)
}

fn read_fasta_bytes(path: &Path) -> Result<Vec<u8>, String> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| format!("{}: {e}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        return Err(format!(
            "{} is gzipped; decompress it first (the corpus is read as plain FASTA)",
            path.display()
        ));
    }
    Ok(raw)
}

/// Which FASTAs the manifest already covers, and how many rows the matrix holds.
fn manifest_state(manifest: &Path) -> (HashSet<String>, usize) {
    let Ok(file) = File::open(manifest) else {
        return (HashSet::new(), 0);
    };
    let (mut done, mut rows) = (HashSet::new(), 0usize);
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if let Some(rest) = line.strip_prefix("#done\t") {
            done.insert(rest.to_string());
        } else if !line.starts_with('#') && !line.is_empty() {
            rows += 1;
        }
    }
    (done, rows)
}

// Hash function to detect duplicates, FNV-1a
fn hash_seq(seq: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in seq {
        h ^= b as u64;
        h = h.wrapping_mul(0x1000_0000_01b3);
    }
    h
}

/// Embed every record of every listed FASTA, one file at a time. Rows go to a raw f32 file,
/// not a `.npy`: that header fixes the row count up front and cannot survive an interrupt.
async fn embed_corpus(args: &EmbedArgs, embedder: &Embedder) -> Result<(PathBuf, usize), String> {
    fs::create_dir_all(&args.work_dir).map_err(|e| e.to_string())?;
    let matrix_path = args.work_dir.join("embeddings.f32");
    let manifest_path = args.work_dir.join("embeddings.ids.tsv");

    let (already_done, mut rows) = if args.resume {
        manifest_state(&manifest_path)
    } else {
        let _ = fs::remove_file(&matrix_path);
        let _ = fs::remove_file(&manifest_path);
        (HashSet::new(), 0)
    };

    let corpus = read_corpus(&args.fasta_list)?;
    let open_append = |p: &Path| {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(p)
            .map_err(|e| format!("{}: {e}", p.display()))
    };
    let mut matrix = BufWriter::new(open_append(&matrix_path)?);
    let mut manifest = BufWriter::new(open_append(&manifest_path)?);
    let mut seen: HashSet<u64> = HashSet::new();
    let mut duplicates = 0usize;

    for (i, entry) in corpus.iter().enumerate() {
        let key = entry.path.display().to_string();
        let tag = format!("[{}/{}]", i + 1, corpus.len());
        if already_done.contains(&key) {
            println!("{tag} {key} — already done");
            continue;
        }

        // This removes duplicates
        let mut records = fasta::parse(&read_fasta_bytes(&entry.path)?)?;
        if !args.keep_duplicates {
            let before = records.len();
            records.retain(|r| seen.insert(hash_seq(&r.seq)));
            duplicates += before - records.len();
        }
        if records.is_empty() {
            println!("{tag} {key} — nothing new");
            writeln!(manifest, "#done\t{key}").map_err(|e| e.to_string())?;
            continue;
        }

        let started = Instant::now();
        // One record is one protein is one embedding row; nothing is pooled per file.
        let (vectors, meta) = embedder
            .embed_records(&records, args.max_len, BatchPolicy::Auto, |_, _| {})
            .await;

        for v in &vectors {
            matrix.write_all(&v.to_le_bytes()).map_err(|e| e.to_string())?;
        }
        for id in &meta.ids {
            writeln!(manifest, "{id}\t{key}\t{}", entry.label).map_err(|e| e.to_string())?;
        }
        writeln!(manifest, "#done\t{key}").map_err(|e| e.to_string())?;
        matrix.flush().map_err(|e| e.to_string())?;
        manifest.flush().map_err(|e| e.to_string())?;

        rows += meta.n;
        println!(
            "{tag} {key} — {} proteins in {:.1} s ({rows} total)",
            meta.n,
            started.elapsed().as_secs_f32(),
        );
    }

    if duplicates > 0 {
        println!("dropped {duplicates} exact-duplicate sequences");
    }
    Ok((matrix_path, rows))
}

/// Read a row-major `[rows, cols]` f32 matrix, checking the length is what it claims.
fn load_f32(path: &Path, rows: usize, cols: usize) -> Result<Vec<f32>, String> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| format!("{}: {e}", path.display()))?;
    let expected = rows * cols * 4;
    if raw.len() != expected {
        return Err(format!(
            "{} is {} bytes, expected {expected} for {rows} x {cols} f32",
            path.display(),
            raw.len(),
        ));
    }
    Ok(raw
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect())
}

fn manifest_rows(manifest: &Path) -> usize {
    manifest_state(manifest).1
}

/// MSE of a trained model over `x`/`y`. f64 accumulation: 800k+ squared terms drift in f32.
fn mse_over<B: Backend>(model: &UmapMlp<B>, x: &[f32], y: &[f32], device: &B::Device) -> f32 {
    let pred = model.project(x, device);
    let sse: f64 = pred
        .iter()
        .zip(y)
        .map(|(p, t)| ((p - t) as f64).powi(2))
        .sum();
    (sse / y.len() as f64) as f32
}

/// Deterministic shuffle, so a rerun holds out the same rows.
fn shuffled_indices(n: usize, seed: u64) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..n).collect();
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    for i in (1..n).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        idx.swap(i, (state % (i as u64 + 1)) as usize);
    }
    idx
}

fn gather(src: &[f32], rows: &[usize], cols: usize) -> Vec<f32> {
    let mut out = Vec::with_capacity(rows.len() * cols);
    for &r in rows {
        out.extend_from_slice(&src[r * cols..(r + 1) * cols]);
    }
    out
}

type Fit = Autodiff<esm_bridge::Gpu>;
type FitCpu = Autodiff<esm_bridge::Cpu>;

/// Largest held-out sample the O(n^2) metric is run on. The estimate is stable long before
/// the whole set is used, and 48k rows would take hours.
const NEIGHBOURHOOD_SAMPLE_MAX: usize = 4_000;

/// Mean fraction of each row's `k` nearest neighbours in 320-D still among its `k` nearest
/// in 2-D — the number that says whether the encoder reproduced the layout. Reads the first
/// `n` rows of `x` and `y`.
fn neighbourhood_preservation(x: &[f32], y: &[f32], n: usize, k: usize) -> f32 {
    let knn = |data: &[f32], dim: usize, cosine: bool| -> Vec<Vec<usize>> {
        let norm: Vec<f32> = if cosine {
            let mut v = data.to_vec();
            for i in 0..n {
                let row = &mut v[i * dim..(i + 1) * dim];
                let m = row.iter().map(|a| a * a).sum::<f32>().sqrt();
                if m > 0.0 {
                    row.iter_mut().for_each(|a| *a /= m);
                }
            }
            v
        } else {
            data.to_vec()
        };
        (0..n)
            .map(|i| {
                let mut d: Vec<(f32, usize)> = (0..n)
                    .filter(|&j| j != i)
                    .map(|j| {
                        let (a, b) = (&norm[i * dim..(i + 1) * dim], &norm[j * dim..(j + 1) * dim]);
                        let v = if cosine {
                            1.0 - a.iter().zip(b).map(|(p, q)| p * q).sum::<f32>()
                        } else {
                            a.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum::<f32>()
                        };
                        (v, j)
                    })
                    .collect();
                d.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                d.into_iter().take(k).map(|(_, j)| j).collect()
            })
            .collect()
    };

    let high = knn(x, HIDDEN, true);
    let low = knn(y, N_DIM, false);
    let total: usize = high
        .iter()
        .zip(&low)
        .map(|(h, l)| {
            let set: HashSet<usize> = h.iter().copied().collect();
            l.iter().filter(|j| set.contains(j)).count()
        })
        .sum();
    total as f32 / (n * k) as f32
}

pub fn main() {
    let result = match Args::parse().command {
        Command::Embed(args) => embed(args),
        Command::Fit(args) => fit(args),
        Command::Eval(args) => eval_encoder(args),
    };
    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn embed(args: EmbedArgs) -> Result<(), String> {
    let weights = fs::read(&args.weights)
        .map_err(|e| format!("could not read {}: {e}", args.weights.display()))?;
    let embedder = if args.cpu {
        Embedder::new_cpu(weights)
    } else {
        pollster::block_on(Embedder::new_gpu(weights, None, false))
    };

    let (_, rows) = pollster::block_on(embed_corpus(&args, &embedder))?;
    println!("corpus: {rows} proteins x {HIDDEN} dimensions");
    println!("next: scripts/umap_coords.py --work-dir {}", args.work_dir.display());
    Ok(())
}

/// Score an encoder that already exists. The split is reproduced from `seed` and `holdout`,
/// so it must match the run that trained it.
fn eval_encoder(args: EvalArgs) -> Result<(), String> {
    let n = manifest_rows(&args.work_dir.join("embeddings.ids.tsv"));
    let x = load_f32(&args.work_dir.join("embeddings.f32"), n, HIDDEN)?;
    let y = load_f32(&args.work_dir.join("coords.f32"), n, N_DIM)?;

    let order = shuffled_indices(n, args.seed);
    let n_hold = ((n as f32 * args.holdout) as usize).clamp(1, n / 2);
    let (hold, train) = order.split_at(n_hold);

    let bytes = fs::read(&args.encoder).map_err(|e| format!("{}: {e}", args.encoder.display()))?;
    let model = UmapMlp::<esm_bridge::Cpu>::from_bytes(bytes, &Default::default())?;
    println!(
        "{}: {} train / {} test rows at holdout {}",
        args.encoder.display(),
        train.len(),
        hold.len(),
        args.holdout,
    );

    let (xt, yt) = (gather(x.as_slice(), train, HIDDEN), gather(&y, train, N_DIM));
    let (xh, yh) = (gather(x.as_slice(), hold, HIDDEN), gather(&y, hold, N_DIM));
    let train_mse = mse_over(&model, &xt, &yt, &Default::default());
    let mse = mse_over(&model, &xh, &yh, &Default::default());
    println!(
        "train MSE {train_mse:.4} vs held-out MSE {mse:.4} — ratio {:.2}",
        mse / train_mse
    );

    let mean = yh.iter().sum::<f32>() / yh.len() as f32;
    let var = yh.iter().map(|t| (t - mean) * (t - mean)).sum::<f32>() / yh.len() as f32;
    println!("held-out MSE {mse:.4} against coordinate variance {var:.4} (R^2 {:.3})", 1.0 - mse / var);

    let n_metric = hold.len().min(NEIGHBOURHOOD_SAMPLE_MAX);
    let k = 15.min((n_metric - 1) / 4).max(1);
    let pred = model.project(&xh[..n_metric * HIDDEN], &Default::default());
    let kept = neighbourhood_preservation(&xh[..n_metric * HIDDEN], &pred, n_metric, k);
    println!(
        "neighbourhood preservation at k={k}: {:.1}% (random would be ~{:.1}%), measured on \
         {n_metric} of {} held-out proteins",
        kept * 100.0,
        100.0 * k as f32 / (n_metric - 1) as f32,
        hold.len(),
    );
    Ok(())
}

/// Fit the encoder by ordinary supervised regression onto the UMAP coordinates.
fn fit(args: FitArgs) -> Result<(), String> {
    let manifest = args.work_dir.join("embeddings.ids.tsv");
    let n = manifest_rows(&manifest);
    if n < 10 {
        return Err(format!("{} lists {n} proteins; too few to fit", manifest.display()));
    }
    let x = load_f32(&args.work_dir.join("embeddings.f32"), n, HIDDEN)?;
    let coords_path = args.work_dir.join("coords.f32");
    if !coords_path.exists() {
        return Err(format!(
            "{} not found; run scripts/umap_coords.py --work-dir {} first",
            coords_path.display(),
            args.work_dir.display(),
        ));
    }
    let y = load_f32(&coords_path, n, N_DIM)?;

    let order = shuffled_indices(n, args.seed);
    let n_hold = ((n as f32 * args.holdout) as usize).clamp(1, n / 2);
    let (hold, train) = order.split_at(n_hold);
    println!("fitting on {} proteins, holding out {}", train.len(), hold.len());

    if args.cpu {
        run_fit::<FitCpu>(&args, &x, &y, train, hold, &Default::default())
    } else {
        run_fit::<Fit>(&args, &x, &y, train, hold, &Default::default())
    }
}

fn run_fit<B: burn::tensor::backend::AutodiffBackend>(
    args: &FitArgs,
    x: &[f32],
    y: &[f32],
    train: &[usize],
    hold: &[usize],
    device: &B::Device,
) -> Result<(), String> {
    let (xt, yt) = (gather(x, train, HIDDEN), gather(y, train, N_DIM));
    let mut model = UmapMlp::<B>::new(device);
    let mut optim = AdamConfig::new().init();
    let started = Instant::now();

    for epoch in 0..args.epochs {
        // Reshuffle each epoch, seeded on the epoch so the run stays reproducible.
        let perm = shuffled_indices(train.len(), args.seed ^ epoch as u64);
        let mut total = 0.0f32;
        let mut steps = 0usize;

        for chunk in perm.chunks(args.batch_size) {
            let xb = Tensor::<B, 2>::from_data(
                TensorData::new(gather(&xt, chunk, HIDDEN), [chunk.len(), HIDDEN]),
                device,
            );
            let yb = Tensor::<B, 2>::from_data(
                TensorData::new(gather(&yt, chunk, N_DIM), [chunk.len(), N_DIM]),
                device,
            );
            let loss = MseLoss::new().forward(model.forward(xb), yb, Reduction::Mean);
            total += loss.clone().into_scalar().to_f32();
            steps += 1;
            let grads = GradientsParams::from_grads(loss.backward(), &model);
            model = optim.step(args.lr, model, grads);
        }
        if epoch % 25 == 0 || epoch + 1 == args.epochs {
            println!("epoch {}/{}: loss {:.6}", epoch + 1, args.epochs, total / steps as f32);
        }
    }
    println!("fitted in {:.1} s", started.elapsed().as_secs_f32());

    // Held-out report. `.valid()` drops the autodiff wrapper for inference and for saving.
    let eval = model.valid();
    let (xh, yh) = (gather(x, hold, HIDDEN), gather(y, hold, N_DIM));
    let pred = eval.project(&xh, &Default::default());

    // Both from the same finished model: the per-epoch loss is a running mean taken while the
    // weights were still moving, so it is not the trained model's error on the training set.
    let train_mse = mse_over(&eval, &xt, &yt, &Default::default());
    let mse = mse_over(&eval, &xh, &yh, &Default::default());
    println!(
        "train MSE {train_mse:.4} vs held-out MSE {mse:.4} — ratio {:.2} ({} train / {} test rows)",
        mse / train_mse,
        yt.len() / N_DIM,
        yh.len() / N_DIM,
    );

    let mean = yh.iter().sum::<f32>() / yh.len() as f32;
    let var = yh.iter().map(|t| (t - mean) * (t - mean)).sum::<f32>() / yh.len() as f32;
    println!("held-out MSE {mse:.4} against coordinate variance {var:.4} (R^2 {:.3})", 1.0 - mse / var);

    // The metric is O(n^2), so it runs on a bounded prefix of the held-out rows, which are
    // already shuffled. At corpus scale the full set would cost hours for a diagnostic.
    let n_metric = hold.len().min(NEIGHBOURHOOD_SAMPLE_MAX);
    // k must stay well under the sample or the measure is trivially 100%. The random
    // baseline is printed so a degenerate run is visible.
    let k = 15.min((n_metric - 1) / 4).max(1);
    // Slice rather than pass the whole held-out set: the metric copies its input.
    let kept = neighbourhood_preservation(
        &xh[..n_metric * HIDDEN],
        &pred[..n_metric * N_DIM],
        n_metric,
        k,
    );
    println!("neighbourhood preservation at k={k}: {:.1}% (random would be ~{:.1}%), \
              measured on {n_metric} of {} held-out proteins",
             kept * 100.0, 100.0 * k as f32 / (n_metric - 1) as f32, hold.len());

    if let Some(parent) = args.out.parent() {
        let _ = fs::create_dir_all(parent);
    }
    // Overwrite: retraining with different parameters must replace the old encoder.
    let mut store = BurnpackStore::from_file(&args.out).overwrite(true);
    eval.save_into(&mut store)
        .map_err(|e| format!("could not write {}: {e}", args.out.display()))?;

    // Read back through the path the browser uses. Tolerance, not equality: the fit may have
    // run on the GPU and this reload is always CPU.
    let written = fs::read(&args.out).map_err(|e| format!("{}: {e}", args.out.display()))?;
    let reloaded = UmapMlp::<esm_bridge::Cpu>::from_bytes(written, &Default::default())?;
    let check = reloaded.project(&xh, &Default::default());
    let worst = check
        .iter()
        .zip(&pred)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    if !(worst < 1e-3) {
        return Err(format!(
            "the saved encoder does not reproduce the trained one: worst coordinate differs by {worst}"
        ));
    }
    println!("wrote {} ({} bytes, reload verified)", args.out.display(),
             fs::metadata(&args.out).map(|m| m.len()).unwrap_or(0));
    Ok(())
}
}

#[cfg(not(target_family = "wasm"))]
fn main() {
    native::main();
}

#[cfg(target_family = "wasm")]
fn main() {}
