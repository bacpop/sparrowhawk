use amr_bridge::{load_index, KmerAssignment};
use std::env;
use std::path::Path;

fn main() -> Result<(), String> {
    let path = env::args()
        .nth(1)
        .ok_or_else(|| "usage: cargo run --example index_stats -- <index.amridx>".to_string())?;
    let index = load_index(Path::new(&path))?;

    let mut gene = 0usize;
    let mut family = 0usize;
    for assignment in index.kmer_map.values() {
        match assignment {
            KmerAssignment::Gene(_) => gene += 1,
            KmerAssignment::Family(_) => family += 1,
        }
    }

    println!(
        "path={path}\nk={}\ngenes={}\nfamilies={}\ntotal_kmers={}\ngene_kmers={}\nfamily_kmers={}",
        index.k,
        index.genes.len(),
        index.families.len(),
        index.kmer_map.len(),
        gene,
        family
    );
    Ok(())
}
