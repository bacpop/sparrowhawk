//! Rewrite the weights burnpack with f16 tensors. Run by `stage_weights.sh`.
//!
//! Lossless here: ESM-2's published weights are already exact f16 values (checked)
//! This is because the training uses mixed-precision weights, and because of that
//! this situation can happen, or if not values of 32bits but with lots of zeros


#[cfg(not(target_family = "wasm"))]
mod native {
    use std::path::PathBuf;

    use burn::tensor::Bytes;
    use burn_store::{BurnpackStore, HalfPrecisionAdapter, ModuleSnapshot};
    use esm_bridge::{model::Model, BatchPolicy, Cpu, CpuDevice, Embedder, HIDDEN};

    const FASTA: &str = "tests/data/proteins20.faa";

    pub fn main() -> Result<(), String> {
        let mut args = std::env::args().skip(1);
        let src: PathBuf = args
            .next()
            .unwrap_or_else(|| "model/generated/esm2_t6_8M_UR50D.bpk".into())
            .into();
        let dst: PathBuf = args
            .next()
            .unwrap_or_else(|| "model/generated/esm2_t6_8M_UR50D.f16.bpk".into())
            .into();

        let raw = std::fs::read(&src).map_err(|e| format!("{}: {e}", src.display()))?;
        let device = CpuDevice::default();

        // Save with the adapter: it reads the snapshot dtype and casts F32 -> F16.
        // And using normal/official burn types!
        let model = Model::<Cpu>::from_bytes(Bytes::from_bytes_vec(raw.clone()), &device);
        let mut out = BurnpackStore::from_file(&dst)
            .overwrite(true)
            .with_to_adapter(HalfPrecisionAdapter::new());
        model
            .save_into(&mut out)
            .map_err(|e| format!("writing {}: {e}", dst.display()))?;

        let halved = std::fs::read(&dst).map_err(|e| format!("{}: {e}", dst.display()))?;
        println!(
            "{} -> {}\n  {} bytes -> {} bytes ({:.0}%)",
            src.display(),
            dst.display(),
            raw.len(),
            halved.len(),
            100.0 * halved.len() as f64 / raw.len() as f64,
        );


        // Now, check
        // Load it back the way the browser will, and require identical output. The same
        // adapter widens F16 -> F32 here; direction is inferred from the stored dtype.
        let records = esm_bridge::fasta::parse(
            &std::fs::read(FASTA).map_err(|e| format!("{FASTA}: {e}"))?,
        )
        .map_err(|e| e.to_string())?;

        // Both go through the ordinary constructor: it detects F16 storage and widens.
        let before = Embedder::new_cpu(raw);
        let after = Embedder::new_cpu(halved);
        let (a, meta) = pollster::block_on(before.embed_records(
            &records,
            1022,
            BatchPolicy::Auto,
            |_, _| {},
        ));
        let (b, _) = pollster::block_on(after.embed_records(
            &records,
            1022,
            BatchPolicy::Auto,
            |_, _| {},
        ));

        let worst = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0f32, f32::max);
        println!(
            "  {} proteins x {HIDDEN} dims, worst coordinate difference {worst:.3e}",
            meta.n
        );
        if worst != 0.0 {
            return Err(format!(
                "f16 storage is not lossless for these weights: worst difference {worst:.3e}. \
                 Ship the f32 burnpack instead."
            ));
        }
        println!("  identical — f16 storage is lossless for these weights");
        Ok(())
    }
}

#[cfg(not(target_family = "wasm"))]
fn main() {
    if let Err(e) = native::main() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

#[cfg(target_family = "wasm")]
fn main() {}
