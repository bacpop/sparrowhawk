//! Dump raw `last_hidden_state` for one sequence, helps for debugging

// Native-only; see the note in smoke.rs.
#[cfg(not(target_family = "wasm"))]
mod native {

use std::fs;
use std::io::Write;

use esm_bridge::{fasta, Embedder, HIDDEN};

pub fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    let fasta_path = args.get(1).map(String::as_str).unwrap_or("tests/data/seqs.faa");
    let index: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
    let out_path = args.get(3).map(String::as_str).unwrap_or("/tmp/hidden_burn.npy");

    let weights = fs::read("model/generated/esm2_t6_8M_UR50D.bpk").map_err(|e| e.to_string())?;
    let records = fasta::parse(&fs::read(fasta_path).map_err(|e| e.to_string())?)?;
    let rec = &records[index];
    eprintln!("sequence: {} ({} raw bytes)", rec.id, rec.seq.len());

    let embedder = Embedder::new_cpu(weights);
    let (hidden, b, l) =
        pollster::block_on(embedder.hidden_states(&[rec.seq.as_slice()], 1022));
    eprintln!("hidden: [{b}, {l}, {HIDDEN}]");

    let rows = b * l;
    let dict = format!(
        "{{'descr': '<f4', 'fortran_order': False, 'shape': ({rows}, {HIDDEN}), }}"
    );
    let padding = (64 - ((10 + dict.len() + 1) % 64)) % 64;
    let header = format!("{dict}{}\n", " ".repeat(padding));

    let mut f = fs::File::create(out_path).map_err(|e| e.to_string())?;
    f.write_all(b"\x93NUMPY").map_err(|e| e.to_string())?;
    f.write_all(&[1u8, 0u8]).map_err(|e| e.to_string())?;
    f.write_all(&(header.len() as u16).to_le_bytes()).map_err(|e| e.to_string())?;
    f.write_all(header.as_bytes()).map_err(|e| e.to_string())?;
    for v in &hidden {
        f.write_all(&v.to_le_bytes()).map_err(|e| e.to_string())?;
    }
    eprintln!("wrote {out_path}");
    Ok(())
}

} // mod native

#[cfg(not(target_family = "wasm"))]
fn main() -> Result<(), String> {
    native::run()
}

#[cfg(target_family = "wasm")]
fn main() {}
