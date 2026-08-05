//! `wasm-bindgen` surface consumed by `www/src/workers/Embedder.ts`.

use wasm_bindgen::prelude::*;

use crate::{fasta, tokenizer::MAX_RESIDUES, BatchPolicy, Embedder, HIDDEN};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_name = postMessage)]
    fn post_message(data: &JsValue);
}

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

fn js_error(error: impl core::fmt::Display) -> JsValue {
    JsValue::from_str(&error.to_string())
}

/// `{ embedProgress: true, done, total }`: for progresses
fn post_progress(done: usize, total: usize) {
    let o = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&o, &"embedProgress".into(), &JsValue::TRUE);
    let _ = js_sys::Reflect::set(&o, &"done".into(), &JsValue::from_f64(done as f64));
    let _ = js_sys::Reflect::set(&o, &"total".into(), &JsValue::from_f64(total as f64));
    post_message(&o.into());
}

/// `{ gpuEvent: true, message }`: for gpu failures
pub fn report_event(message: &str) {
    let o = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&o, &"gpuEvent".into(), &JsValue::TRUE);
    let _ = js_sys::Reflect::set(&o, &"message".into(), &JsValue::from_str(message));
    post_message(&o.into());
}

/// One FASTA's results, as `Float32Array`s so the worker can transfer rather than clone.
#[wasm_bindgen]
pub struct EmbeddingBatch {
    meta: String,
    vectors: Vec<f32>,
    coords: Vec<f32>,
}

#[wasm_bindgen]
impl EmbeddingBatch {
    /// JSON: `{ ids, desc, lengths, truncated, dim, n, backend, pooling }`.
    #[wasm_bindgen(getter)]
    pub fn meta(&self) -> String {
        self.meta.clone()
    }

    /// Row-major `[n, 320]`.
    #[wasm_bindgen(getter)]
    pub fn vectors(&self) -> js_sys::Float32Array {
        js_sys::Float32Array::from(&self.vectors[..])
    }

    /// Row-major `[n, 2]` from the UMAP encoder, or empty if none was supplied.
    #[wasm_bindgen(getter)]
    pub fn coords(&self) -> js_sys::Float32Array {
        js_sys::Float32Array::from(&self.coords[..])
    }
}

/// Required: the coordinates it produces are part of every result, so absence is fatal.
fn with_encoder(mut inner: Embedder, encoder: Option<Vec<u8>>) -> Result<Embedder, JsValue> {
    let bytes = encoder.ok_or_else(|| js_error("no UMAP encoder was supplied"))?;
    inner.set_encoder(bytes).map_err(js_error)?;
    Ok(inner)
}

// Main struct for the wasm/js interface, like in the assembler, in ska, etc.
#[wasm_bindgen]
pub struct EsmEmbedder {
    inner: Embedder,
}

#[wasm_bindgen]
impl EsmEmbedder {
    /// Async factory: `wasm-bindgen` constructors cannot be async. `Err` means retry with
    /// `use_gpu = false`.
    pub async fn create(
        weights: Vec<u8>,
        encoder: Option<Vec<u8>>,
        use_gpu: bool,
        power_pref: u32,
        tasks_max: usize,
    ) -> Result<EsmEmbedder, JsValue> {
        init_panic_hook();

        if use_gpu {
            let tasks_max = if tasks_max == 0 { None } else { Some(tasks_max) };
            // 2 is the low-power entry; see `listGpuAdapters` in platform/gpu.ts.
            let inner = Embedder::new_gpu(weights, tasks_max, power_pref == 2).await;
            crate::report_event("model weights uploaded, running self-check");
            // Getting this far proves nothing: WebGPU creates the device, then rejects
            // shaders. Verify before the user's first file.
            inner.warmup_checked().await.map_err(js_error)?;
            crate::report_event("self-check passed, GPU is computing correctly");
            return Ok(EsmEmbedder { inner: with_encoder(inner, encoder)? });
        }

        Ok(EsmEmbedder {
            inner: with_encoder(Embedder::new_cpu(weights), encoder)?,
        })
    }

    /// JSON description of the loaded model, for the UI.
    pub fn info(&self) -> String {
        self.inner.info()
    }

    /// `"webgpu"` or `"cpu"`.
    pub fn backend(&self) -> String {
        self.inner.backend().to_string()
    }

    /// Embedding dimensionality (320 for esm2_t6_8M_UR50D).
    pub fn dim(&self) -> usize {
        HIDDEN
    }

    /// Embed a protein FASTA; `fasta_bytes` must already be decompressed. Proteins longer
    /// than [`MAX_RESIDUES`] are head-truncated and flagged in the metadata.
    pub async fn embed_fasta(&self, fasta_bytes: Vec<u8>) -> Result<EmbeddingBatch, JsValue> {
        let records = fasta::parse(&fasta_bytes).map_err(js_error)?;
        let (vectors, meta) = self
            .inner
            .embed_records(&records, MAX_RESIDUES, BatchPolicy::Auto, post_progress)
            .await;
        let coords = self.inner.project(&vectors);
        let meta = serde_json::to_string(&meta).map_err(js_error)?;
        Ok(EmbeddingBatch { meta, vectors, coords })
    }
}
