//! ESM-2 (8M) protein embeddings.

// burn-onnx emits no_std-friendly code that refers to `alloc::` paths directly.
extern crate alloc;

pub mod fasta;
pub mod projector;
pub mod tokenizer;

/// Generated at build time by burn-onnx from `model/esm2_t6_8M_UR50D.onnx`.
#[allow(clippy::all, dead_code, unused_imports)]
pub mod model {
    include!(concat!(env!("OUT_DIR"), "/model/esm2_t6_8M_UR50D.rs"));
}

use burn::prelude::*;
use burn::tensor::{Bytes, Int, Tensor, TensorData};
use burn_store::{BurnpackStore, HalfPrecisionAdapter, ModuleSnapshot};
use model::Model;
use tokenizer::MAX_RESIDUES;

/// Hidden/embedding size of esm2_t6_8M_UR50D.
pub const HIDDEN: usize = 320;
const N_HEADS: usize = 20; // Attention heads of ESM

// Attention scores for one batch, in f32 elements (x4 for bytes). Only reached on a device
// reporting near WebGPU's guaranteed 128 MiB binding; real adapters report far more.
const ATTN_ELEM_BUDGET_FLOOR: usize = 16_000_000; // i.e. 16 * 10^6 f32 elements, in bytes

// 256 MB. This, not the adapter, is what caps batches in practice: a 2 GiB binding limit
// derives 201M elements. The old 32M was set where wider batches were free; here they are not.
const ATTN_ELEM_BUDGET_CEILING: usize = 64_000_000; // i.e. 32 * 10^6 f32 elements, in bytes

// NdArray allocates on the heap, so address space is the limit, not a binding size. The
// browser gets less: a worker sharing a 4 GiB wasm heap is not a desktop.
#[cfg(target_family = "wasm")]
const ATTN_ELEM_BUDGET_CPU: usize = 32_000_000;
#[cfg(not(target_family = "wasm"))]
const ATTN_ELEM_BUDGET_CPU: usize = 64_000_000;

// Dispatches per queue.submit, when the caller does not choose. 1, not cubecl's 32: the i915
// timeout is per *submission*, although at the end it seems it has been fixed (?)
const GPU_TASKS_MAX: usize = 1;

// Only the 128-token bin reaches this; above it the budget (of how many elements) binds first.
const MAX_BATCH_GPU: usize = 64;

/// NdArray gains little from wide batches while peak memory rises.
const MAX_BATCH_CPU: usize = 8;

/// Probe sequence for [`Embedder::warmup_checked`].
const WARMUP_SEQ: &[u8] = b"MK";

/// Known-good first elements of the embedding of [`WARMUP_SEQ`], pinned to onnxruntime by `tests/fidelity.rs`.
const WARMUP_EXPECTED: [f32; 8] = [
    0.076_747_07,
    -0.184_892_62,
    0.428_570_5,
    0.382_838_2,
    -0.132_557_48,
    -0.145_408_18,
    -0.344_398_62,
    -0.148_791_76,
];

/// Comfortably above GPU float-reordering noise (~4e-5), well below a real failure.
const WARMUP_TOLERANCE: f32 = 1e-2;

use core::sync::atomic::{AtomicBool, Ordering};

// cubecl keeps one global client per device; init_setup_async panics on a second call.
static GPU_INITIALISED: AtomicBool = AtomicBool::new(false);

pub type Cpu = burn::backend::NdArray;
pub type CpuDevice = burn::backend::ndarray::NdArrayDevice;
pub type Gpu = burn::backend::WebGpu;
pub type GpuDevice = burn::backend::wgpu::WgpuDevice;

/// Production is always `Auto`; `Fixed` is for tests and benchmark sweeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchPolicy {
    Auto,
    Fixed(usize),
}

/// Per-sequence metadata accompanying the embedding matrix.
#[derive(Debug, Clone, serde::Serialize)]
pub struct EmbeddingMeta {
    pub ids: Vec<String>,
    pub desc: Vec<String>,
    /// Residues actually embedded (post-truncation).
    pub lengths: Vec<usize>,
    pub truncated: Vec<bool>,
    pub dim: usize,
    pub n: usize,
    pub backend: String,
    pub pooling: &'static str,

    /// Smallest and largest batch the planner used, so the automatic choice is visible.
    pub batch_min: usize,
    pub batch_max: usize,

    /// Attention-score budget in f32 elements, resolved from the device.
    pub budget_elems: usize,
}

/// The generic inference kernel: row-major `[batch, HIDDEN]`, mask-aware mean-pooled.
async fn forward_batch<B: Backend>(
    model: &Model<B>,
    device: &B::Device,
    enc: &tokenizer::Encoded,
) -> Vec<f32> {
    let (b, l) = (enc.batch, enc.len);

    // reshapes are not needed, we could remove them, but they are apparently inoquous as dimensions changes are checked when on-device
    let ids = Tensor::<B, 1, Int>::from_data(TensorData::new(enc.input_ids.clone(), [b * l]), device)
        .reshape([b, l]);

    let att =
        Tensor::<B, 1, Int>::from_data(TensorData::new(enc.attention_mask.clone(), [b * l]), device)
            .reshape([b, l]);


    // hidden is a tensor with size [batch size, length (padded) of this batch, embedding vector]
    let hidden = model.forward(ids, att);

    // Mean over residues only: pool_mask is zero on <cls>, <eos> and <pad>.
    let pool = Tensor::<B, 1>::from_data(TensorData::new(enc.pool_mask.clone(), [b * l]), device)
        .reshape([b, l]);

    // clamp_min guards the degenerate all-padding row (an empty FASTA record).
    let denom = pool.clone().sum_dim(1).clamp_min(1.0);

    let pooled = (hidden * pool.unsqueeze_dim::<3>(2))
        .sum_dim(1)
        .squeeze_dim::<2>(1)
        / denom;

    pooled
        .into_data_async()
        .await
        .expect("reading pooled embeddings back from the device failed")
        .convert::<f32>()
        .into_vec::<f32>()
        .expect("pooled embeddings should convert to Vec<f32>")
}

/// Without these every failure mode looks identical: `BufferAsyncError` is an empty struct
/// on the web (wgpu #2939).
fn install_gpu_diagnostics(setup: &burn::backend::wgpu::WgpuSetup) {
    setup.device.set_device_lost_callback(|reason, message| {
        report_event(&format!("device lost ({reason:?}): {message}"));
    });
    setup
        .device
        .on_uncaptured_error(std::sync::Arc::new(|error| {
            report_event(&format!("uncaptured GPU error: {error}"));
        }));
}

/// The worker's message channel in the browser, stderr natively. Not GPU-only: the batch log
/// below runs on either backend, so the prefix must not claim otherwise.
pub(crate) fn report_event(message: &str) {
    #[cfg(target_family = "wasm")]
    crate::wasm::report_event(message);
    #[cfg(not(target_family = "wasm"))]
    eprintln!("[esm-bridge] {message}");
}

/// Read from the adapter, not assumed. The headroom is deliberate: `BufferTooBig` panics.
fn gpu_budget_elems(device: &GpuDevice) -> usize {
    use burn::backend::wgpu::WgpuRuntime;
    use burn::cubecl::Runtime;

    let client = WgpuRuntime::client(device);
    let max_bytes = client.properties().memory.max_page_size;
    let usable = (max_bytes / 4) as usize / 2 * 3 / 4; // elements, halved, then 3/4
    let budget = usable.clamp(ATTN_ELEM_BUDGET_FLOOR, ATTN_ELEM_BUDGET_CEILING);

    // The peak figure is irreducible: one max-length sequence needs it whatever the batch.
    let peak_mb = N_HEADS * 1024 * 1024 * 4 / (1024 * 1024);
    // This is the number of attention heads times the entries of the attention matrix times 4 because 32 bits are four bytes, and over 1024^2 because of the conversion to MB, or MiB: I always confuse them.
    report_event(&format!(
        "device max binding {} MB, attention budget {} M elements ({} MB), \
         peak single-sequence buffer {peak_mb} MB at {MAX_RESIDUES} residues",
        max_bytes / (1024 * 1024),
        budget / 1_000_000,
        budget * 4 / (1024 * 1024),
    ));
    budget
}

/// `wasm-bindgen` cannot export generics, so the backend is picked here instead.
enum Runner {
    Cpu(Model<Cpu>, CpuDevice),
    Gpu(Model<Gpu>, GpuDevice),
}

/// Dispatch to the active backend, monomorphising `$body` once per arm.
macro_rules! with_runner {
    ($me:expr, |$model:ident, $device:ident| $body:block) => {
        match &$me.runner {
            Runner::Cpu($model, $device) => $body,
            Runner::Gpu($model, $device) => $body,
        }
    };
}

/// Backend-agnostic embedder.
pub struct Embedder {
    runner: Runner,
    backend: &'static str,
    /// Attention-score budget in f32 elements, resolved once at construction.
    budget_elems: usize,
    /// Ceiling the automatic planner will not exceed.
    max_batch: usize,
    /// Always CPU: a 3-layer MLP over a few thousand rows is microseconds on NdArray.
    projector: Option<projector::UmapMlp<Cpu>>,
}

/// True if the burnpack stores any tensor as F16 (`0x63 F 1 6`, a CBOR text string). Checked,
/// not assumed: [`HalfPrecisionAdapter`] reads the stored dtype and would cast an F32 file down.
fn weights_are_half(bytes: &[u8]) -> bool {
    const MARKER: [u8; 4] = [0x63, b'F', b'1', b'6'];
    // The header is tiny next to the payload; bound the scan so this stays O(header).
    let end = bytes.len().min(1 << 20);
    bytes[..end].windows(MARKER.len()).any(|w| w == MARKER)
}

/// Load the generated model, widening F16 storage back to F32 parameters if needed.
fn load_model<B: Backend>(weights: Vec<u8>, device: &B::Device) -> Model<B> {
    let half = weights_are_half(&weights);
    let mut model = Model::<B>::new(device);
    let mut store = BurnpackStore::from_bytes(Some(Bytes::from_bytes_vec(weights)));
    if half {
        store = store.with_from_adapter(HalfPrecisionAdapter::new());
    }
    model
        .load_from(&mut store)
        .expect("failed to load the ESM-2 burnpack");
    model
}

impl Embedder {
    /// Constructor for CPU (NdArray) backend. `weights` is the raw, already-decompressed burnpack.
    pub fn new_cpu(weights: Vec<u8>) -> Self {
        let device = CpuDevice::default();
        let model = load_model::<Cpu>(weights, &device);
        Self {
            runner: Runner::Cpu(model, device),
            backend: "cpu",
            budget_elems: ATTN_ELEM_BUDGET_CPU,
            max_batch: MAX_BATCH_CPU,
            projector: None,
        }
    }

    /// Attach the trained UMAP encoder. Required: callers treat failure as fatal.
    pub fn set_encoder(&mut self, encoder: Vec<u8>) -> Result<(), String> {
        self.projector = Some(projector::UmapMlp::from_bytes(encoder, &CpuDevice::default())?);
        Ok(())
    }

    /// Row-major `[n, 2]` for a row-major `[n, HIDDEN]` matrix, or empty with no encoder.
    pub fn project(&self, vectors: &[f32]) -> Vec<f32> {
        match &self.projector {
            Some(p) => p.project(vectors, &CpuDevice::default()),
            None => Vec::new(),
        }
    }

    /// Attention-score budget in f32 elements, and the batch ceiling.
    pub fn budget_elems(&self) -> usize {
        self.budget_elems
    }
    pub fn max_batch(&self) -> usize {
        self.max_batch
    }

    /// WebGPU backend; any failure means "fall back to [`Embedder::new_cpu`]".
    pub async fn new_gpu(
        weights: Vec<u8>,
        tasks_max: Option<usize>,
        prefer_integrated: bool,
    ) -> Self {
        use burn::backend::wgpu::{init_setup_async, RuntimeOptions, WgpuDevice};
        // Only WebGpu exists in the browser; natively let cubecl pick Vulkan/Metal/DX12.
        #[cfg(target_family = "wasm")]
        use burn::backend::wgpu::graphics::WebGpu as SelectedApi;
        #[cfg(not(target_family = "wasm"))]
        use burn::backend::wgpu::graphics::AutoGraphicsApi as SelectedApi;

        // cubecl maps IntegratedGpu to PowerPreference::LowPower on the web, and to real
        // adapter enumeration natively; anything else asks for high performance.
        let device = if prefer_integrated {
            WgpuDevice::IntegratedGpu(0)
        } else {
            WgpuDevice::DefaultDevice
        };
        // A second init_setup_async panics with "Service already initialized", which the
        // GPU toggle reaches in normal use.
        if GPU_INITIALISED
            .compare_exchange(false, true, Ordering::SeqCst, Ordering::SeqCst)
            .is_ok()
        {
            // A struct literal, not `default()`: that reads an env var, useless in a browser.
            // TODO: adapt for proper native inference
            let tasks_max = tasks_max.unwrap_or(GPU_TASKS_MAX).max(1);
            report_event(&format!("initialising GPU with tasks_max = {tasks_max}"));
            let options = RuntimeOptions {
                tasks_max,
                ..RuntimeOptions::default()
            };
            let setup = init_setup_async::<SelectedApi>(&device, options).await;
            install_gpu_diagnostics(&setup);
        }
        let budget_elems = gpu_budget_elems(&device);
        let model = load_model::<Gpu>(weights, &device);
        Self {
            runner: Runner::Gpu(model, device),
            backend: "webgpu",
            budget_elems,
            max_batch: MAX_BATCH_GPU,
            projector: None,
        }
    }

    pub fn backend(&self) -> &'static str {
        self.backend
    }

    pub fn info(&self) -> String {
        format!(
            r#"{{"model":"facebook/esm2_t6_8M_UR50D","layers":6,"hidden":{HIDDEN},"heads":{N_HEADS},"vocab":{},"max_residues":{MAX_RESIDUES},"backend":"{}"}}"#,
            tokenizer::VOCAB_SIZE,
            self.backend
        )
    }

    /// Embed already-parsed records. `progress` is called with `(done, total)`.
    pub async fn embed_records(
        &self,
        records: &[fasta::Record],
        max_len: usize,
        policy: BatchPolicy,
        mut progress: impl FnMut(usize, usize),
    ) -> (Vec<f32>, EmbeddingMeta) {
        let max_res = max_len.clamp(1, MAX_RESIDUES);
        let n = records.len();

        let counts: Vec<usize> = records
            .iter()
            .map(|r| tokenizer::residue_count(&r.seq).min(max_res))
            .collect();

        // Fixed still needs the memory clamp: an oversized buffer panics rather than erroring.
        let groups = match policy {
            BatchPolicy::Auto => tokenizer::plan_batches_auto(
                &counts,
                self.budget_elems,
                N_HEADS,
                self.max_batch,
            ),
            BatchPolicy::Fixed(size) => tokenizer::plan_batches_fixed(&counts, size.max(1))
                .into_iter()
                .flat_map(|group| {
                    let padded =
                        tokenizer::padded_len(group.iter().map(|&i| counts[i]).max().unwrap_or(0));
                    let cap =
                        tokenizer::batch_cap(padded, self.budget_elems, N_HEADS, self.max_batch);
                    group
                        .chunks(cap)
                        .map(|c| c.to_vec())
                        .collect::<Vec<_>>()
                })
                .collect(),
        };

        let mut vectors = vec![0f32; n * HIDDEN];
        let mut truncated = vec![false; n];
        let mut done = 0usize;
        let mut batch_min = usize::MAX;
        let mut batch_max = 0usize;

        for (gi, group) in groups.iter().enumerate() {
            batch_min = batch_min.min(group.len());
            batch_max = batch_max.max(group.len());

            let seqs: Vec<&[u8]> = group.iter().map(|&i| records[i].seq.as_slice()).collect();
            let enc = tokenizer::encode_batch(&seqs, max_res);

            // Logged before the work, so a GPU death leaves the guilty shape as the last line.
            report_event(&format!(
                "batch {}/{}: {} x {} tokens ({} MB attention)",
                gi + 1,
                groups.len(),
                enc.batch,
                enc.len,
                enc.batch * N_HEADS * enc.len * enc.len * 4 / (1024 * 1024),
            ));

            let out = with_runner!(self, |m, d| { forward_batch(m, d, &enc).await });

            for (j, &idx) in group.iter().enumerate() {
                truncated[idx] = enc.truncated[j];
                vectors[idx * HIDDEN..(idx + 1) * HIDDEN]
                    .copy_from_slice(&out[j * HIDDEN..(j + 1) * HIDDEN]);
            }
            done += group.len();
            progress(done, n);
        }

        let meta = EmbeddingMeta {
            ids: records.iter().map(|r| r.id.clone()).collect(),
            desc: records.iter().map(|r| r.desc.clone()).collect(),
            lengths: counts,
            truncated,
            dim: HIDDEN,
            n,
            backend: self.backend.to_string(),
            pooling: "mean over residues (cls/eos/pad excluded)",
            batch_min: if groups.is_empty() { 0 } else { batch_min },
            batch_max,
            budget_elems: self.budget_elems,
        };
        (vectors, meta)
    }

    /// Parse a FASTA and embed it.
    pub async fn embed_fasta_bytes(
        &self,
        fasta_bytes: &[u8],
        max_len: usize,
        policy: BatchPolicy,
        progress: impl FnMut(usize, usize),
    ) -> Result<(Vec<f32>, EmbeddingMeta), String> {
        let records = fasta::parse(fasta_bytes)?;
        Ok(self.embed_records(&records, max_len, policy, progress).await)
    }

    /// Raw `last_hidden_state`, row-major `[batch * len, HIDDEN]`, for `debug_hidden.rs`.
    pub async fn hidden_states(
        &self,
        seqs: &[&[u8]],
        max_len: usize,
    ) -> (Vec<f32>, usize, usize) {
        let enc = tokenizer::encode_batch(seqs, max_len.clamp(1, MAX_RESIDUES));
        let (b, l) = (enc.batch, enc.len);
        let out = with_runner!(self, |m, d| {
            let ids = Tensor::<_, 1, Int>::from_data(
                TensorData::new(enc.input_ids.clone(), [b * l]),
                d,
            )
            .reshape([b, l]);
            let att = Tensor::<_, 1, Int>::from_data(
                TensorData::new(enc.attention_mask.clone(), [b * l]),
                d,
            )
            .reshape([b, l]);
            m.forward(ids, att)
                .into_data_async()
                .await
                .expect("device readback failed")
                .convert::<f32>()
                .into_vec::<f32>()
                .expect("hidden states should convert to Vec<f32>")
        });
        (out, b, l)
    }

    /// Compile shaders and verify the backend actually computes. WebGPU creates a device
    /// happily and only then rejects shaders, surfacing nothing but wrong numbers.
    pub async fn warmup_checked(&self) -> Result<(), String> {
        // Probe 1: the pinned reference vector, at batch 1 / length 128 after quantisation.
        self.check_probe().await?;

        // Probe 2: probe 1 is 2 residues, so it never compiles the kernels real work uses.
        // Batch 4 against batch 1 exercises those without needing a second pinned constant.
        let long: Vec<u8> = WARMUP_SEQ.iter().copied().cycle().take(300).collect();
        let seq = long.as_slice();
        let batched = tokenizer::encode_batch(&[seq, seq, seq, seq], MAX_RESIDUES);
        let single = tokenizer::encode_batch(&[seq], MAX_RESIDUES);
        let a = with_runner!(self, |m, d| { forward_batch(m, d, &batched).await });
        let b = with_runner!(self, |m, d| { forward_batch(m, d, &single).await });
        if !a.iter().all(|v| v.is_finite()) {
            return Err("backend self-check produced non-finite output at batch 4".into());
        }
        for (i, (&got, &want)) in a.iter().take(HIDDEN).zip(b.iter()).enumerate() {
            if (got - want).abs() > WARMUP_TOLERANCE {
                return Err(format!(
                    "backend self-check failed at dimension {i} for a 300-residue batch of 4: \
                     got {got:.6}, batch-of-1 gives {want:.6}. The backend initialised but is \
                     not computing correctly at realistic shapes."
                ));
            }
        }
        Ok(())
    }

    async fn check_probe(&self) -> Result<(), String> {
        let enc = tokenizer::encode_batch(&[WARMUP_SEQ], MAX_RESIDUES);
        let out = with_runner!(self, |m, d| { forward_batch(m, d, &enc).await });

        if out.len() < WARMUP_EXPECTED.len() {
            return Err(format!(
                "backend self-check returned {} values, expected at least {}",
                out.len(),
                WARMUP_EXPECTED.len()
            ));
        }
        if !out.iter().all(|v| v.is_finite()) {
            return Err("backend self-check produced non-finite output".into());
        }
        for (i, (&got, &want)) in out.iter().zip(WARMUP_EXPECTED.iter()).enumerate() {
            if (got - want).abs() > WARMUP_TOLERANCE {
                return Err(format!(
                    "backend self-check failed at dimension {i}: got {got:.6}, expected \
                     {want:.6}. The backend initialised but is not computing correctly."
                ));
            }
        }
        Ok(())
    }
}

#[cfg(target_family = "wasm")]
pub mod wasm;
