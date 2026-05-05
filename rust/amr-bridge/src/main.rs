use amr_bridge::{
    build_index_from_resfinder_db, compare_to_baseline, detect_fasta, load_index,
    parse_resfinder_json, run_native_resfinder, save_index, DetectParams,
};
use clap::{Parser, Subcommand};
use std::fs;
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    BuildIndex {
        #[arg(long)]
        db_root: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value_t = 31)]
        k: usize,
    },
    Detect {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        fasta: PathBuf,
        #[arg(long)]
        sample_name: Option<String>,
        #[arg(long, default_value_t = 0.05)]
        min_gene_fraction: f64,
        #[arg(long, default_value_t = 0.30)]
        min_family_fraction: f64,
    },
    Eval {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        fasta: PathBuf,
        #[arg(long)]
        resfinder_root: PathBuf,
        #[arg(long)]
        db_root: PathBuf,
        #[arg(long)]
        out_dir: PathBuf,
        #[arg(long, default_value = "/usr/bin/blastn")]
        blastn_path: PathBuf,
        #[arg(long)]
        kma_path: Option<PathBuf>,
        #[arg(long)]
        db_root_kma: Option<PathBuf>,
        #[arg(long)]
        baseline_json: Option<PathBuf>,
    },
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    match cli.command {
        Command::BuildIndex { db_root, out, k } => {
            let index = build_index_from_resfinder_db(&db_root, k)?;
            save_index(&index, &out)?;
            println!(
                "Built index: k={} genes={} families={} kmers={}",
                index.k,
                index.genes.len(),
                index.families.len(),
                index.kmer_map.len()
            );
        }
        Command::Detect {
            index,
            fasta,
            sample_name,
            min_gene_fraction,
            min_family_fraction,
        } => {
            let index = load_index(&index)?;
            let bytes =
                fs::read(&fasta).map_err(|e| format!("Read FASTA {}: {e}", fasta.display()))?;
            let result = detect_fasta(
                &index,
                &bytes,
                sample_name.as_deref().unwrap_or_else(|| {
                    fasta
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample")
                }),
                &DetectParams {
                    min_gene_fraction,
                    min_family_fraction,
                },
            )?;
            println!(
                "{}",
                serde_json::to_string_pretty(&result).map_err(|e| e.to_string())?
            );
        }
        Command::Eval {
            index,
            fasta,
            resfinder_root,
            db_root,
            out_dir,
            blastn_path,
            kma_path,
            db_root_kma,
            baseline_json,
        } => {
            let baseline = if let Some(baseline_json) = baseline_json {
                baseline_json
            } else {
                run_native_resfinder(
                    &resfinder_root,
                    &fasta,
                    &db_root,
                    &out_dir,
                    &blastn_path,
                    kma_path.as_deref(),
                    db_root_kma.as_deref(),
                )?
            };
            let baseline_hits = parse_resfinder_json(&baseline)?;
            let index = load_index(&index)?;
            let bytes =
                fs::read(&fasta).map_err(|e| format!("Read FASTA {}: {e}", fasta.display()))?;
            let sample_name = fasta
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sample");
            let result = detect_fasta(&index, &bytes, sample_name, &DetectParams::default())?;
            let report = compare_to_baseline(&result, &baseline_hits);
            println!(
                "{}",
                serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
            );
        }
    }
    Ok(())
}
