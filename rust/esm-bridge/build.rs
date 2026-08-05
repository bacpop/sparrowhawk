use burn_onnx::{LoadStrategy, ModelGen};
use std::path::PathBuf;

const ONNX: &str = "model/esm2_t6_8M_UR50D.onnx";

fn main() {
    println!("cargo:rerun-if-changed={ONNX}");

    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    if !manifest.join(ONNX).exists() {
        panic!("Missing {ONNX}. Run scripts/fetch_model.sh to obtain it.");
    }

    ModelGen::new()
        .input(ONNX)
        .out_dir("model/")
        // Generate only `Model::from_bytes()`; the weights ship as a separate asset.
        .load_strategy(LoadStrategy::Bytes)
        // MUST stay false: onnx-ir's `coalesce_attention` double-scales ESM-2's pre-scaled
        // Q, squaring the score scale. Upstream tracel-ai/burn-onnx#436, open.
        .simplify(false)
        // Flip to true by hand to emit .onnx.txt / .graph.txt for codegen debugging.
        .development(false)
        .run_from_script();

    let out_dir = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    rewrite_i64_casts(&out_dir.join("model/esm2_t6_8M_UR50D.rs"));

    // Out of the hashed $OUT_DIR so stage_weights.sh can find it.
    let out = out_dir.join("model/esm2_t6_8M_UR50D.bpk");
    if out.exists() {
        let dst = manifest.join("model/generated");
        std::fs::create_dir_all(&dst).ok();
        std::fs::copy(&out, dst.join("esm2_t6_8M_UR50D.bpk")).ok();
    }
}

/// WGSL has no 64-bit integer, so every i64 kernel fails to compile in a browser. This method just changes that to i32. After checking, there was no lost in precision comparing both model's inference
fn rewrite_i64_casts(generated: &std::path::Path) {
    let src = std::fs::read_to_string(generated)
        .unwrap_or_else(|e| panic!("reading generated model {}: {e}", generated.display()));

    let found = src.matches("burn::tensor::DType::I64").count();
    let patched = src.replace("burn::tensor::DType::I64", "burn::tensor::DType::I32");

    // Fail the build rather than ship something no browser can run.
    let remaining = patched.matches("DType::I64").count();
    assert_eq!(
        remaining, 0,
        "generated model still contains {remaining} I64 dtype(s) after rewriting {found}; \
         WebGPU cannot execute those kernels. Inspect {} and extend rewrite_i64_casts.",
        generated.display()
    );

    std::fs::write(generated, patched).expect("writing patched model");
    println!("cargo:warning=esm-bridge: rewrote {found} I64 casts to I32 for WebGPU (WGSL has no i64)");
}
