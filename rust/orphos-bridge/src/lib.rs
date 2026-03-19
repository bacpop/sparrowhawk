use orphos_core::config::{OrphosConfig, OutputFormat};
use orphos_core::engine::{OrphosAnalyzer, UntrainedOrphos};
use orphos_core::sequence::encoded::EncodedSequence;
use orphos_core::types::Training;
use orphos_core::output::write_results;
use orphos_core::results::OrphosResults;
use wasm_bindgen::prelude::*;
use wasm_bindgen_file_reader::WebSysFile;

mod htslib;
use htslib::{bgzf_compress, csi_index_gff, faidx_index_fasta};

pub mod fastx_wasm;
extern crate console_error_panic_hook;

use crate::fastx_wasm::open_fasta;
use flate2::read::MultiGzDecoder;
use seq_io::fasta::Record;
use std::io::{Cursor, Read};

const MIN_NT_CONTIG: usize = 96;

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen]
/// Struct to interact with JS when working with WebAssembly
pub struct OrphosData {
    metag: bool,
    format: OutputFormat,
    closed_ends: bool,
    mask_n_runs: bool,
    force_non_sd: bool,
    translation_table: u8,
    results: Option<Vec<OrphosResults>>,
    gene_count: Option<usize>,
    sequence_count: Option<usize>,
    raw_fasta_bytes: Vec<u8>,
    fasta_bgz: Vec<u8>,
    fasta_fai: Vec<u8>,
    fasta_gzi: Vec<u8>,
    gff_bgz: Vec<u8>,
    gff_csi: Vec<u8>,
}

#[wasm_bindgen]
impl OrphosData {
    /// Constructor
    pub fn new(
        metag: bool,
        format: String,
        closed_ends: bool,
        mask_n_runs: bool,
        force_non_sd: bool,
        translation_table: u8,
    ) -> Self {
        if cfg!(debug_assertions) {
            init_panic_hook();
        }
        wasm_logger::init(wasm_logger::Config::default());

        if !(0..=25).contains(&translation_table)
            || translation_table == 7
            || translation_table == 8
            || (17..=20).contains(&translation_table)
        {
            panic!("Invalid translation table specified");
        }

        let output_format = match format.as_str() {
            "gbk" | "genbank" => OutputFormat::Genbank,
            "gff" => OutputFormat::Gff,
            "sco" => OutputFormat::Sco,
            "gca" => OutputFormat::Gca,
            _ => panic!("Invalid output format"),
        };

        OrphosData {
            metag,
            format: output_format,
            closed_ends,
            mask_n_runs,
            force_non_sd,
            translation_table,
            results: None,
            gene_count: None,
            sequence_count: None,
            raw_fasta_bytes: Vec::new(),
            fasta_bgz: Vec::new(),
            fasta_fai: Vec::new(),
            fasta_gzi: Vec::new(),
            gff_bgz: Vec::new(),
            gff_csi: Vec::new(),
        }
    }

    /// Step 2: read and decompress FASTA file into raw_fasta_bytes
    pub fn read_fasta(&mut self, input_file: web_sys::File) {
        let mut wsf = WebSysFile::new(input_file);
        let mut raw = Vec::new();
        wsf.read_to_end(&mut raw).expect("read failed");
        self.raw_fasta_bytes = if raw.starts_with(&[0x1F, 0x8B]) {
            let mut dec = MultiGzDecoder::new(Cursor::new(&raw));
            let mut out = Vec::new();
            dec.read_to_end(&mut out).expect("gz decode failed");
            out
        } else {
            raw
        };
        log::info!("FASTA read complete: {} bytes", self.raw_fasta_bytes.len());
    }

    /// Step 3: bgzf-compress and faidx-index stored FASTA bytes
    pub fn index_fasta(&mut self) {
        log::info!("Compressing FASTA to BGZF...");
        self.fasta_bgz = bgzf_compress(&self.raw_fasta_bytes);
        log::info!("Indexing FASTA (faidx)...");
        let faidx = faidx_index_fasta(&self.fasta_bgz);
        self.fasta_fai = faidx.fai;
        self.fasta_gzi = faidx.gzi;
        log::info!("FASTA indexing complete.");
    }

    /// Steps 4+5: run orphos; clears raw_fasta_bytes when done
    pub fn call_genes(&mut self) {
        log::info!("Analysing FASTA file. Getting config struct.");
        let orphosconfig = self.orphos_config();
        let mut analyser = OrphosAnalyzer::new(orphosconfig);
        let mut all_results = Vec::new();

        if !self.metag {
            // ── Single mode: train on whole genome, then analyze per-contig ──
            log::info!("Single mode: training pass...");
            let mut training_seq: Vec<u8> = Vec::new();
            let mut cursor = Cursor::new(&self.raw_fasta_bytes);
            let mut reader = open_fasta(&mut cursor);
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid FASTA record");
                let tmpvec = seqrec.full_seq().to_vec();
                if tmpvec.len() >= MIN_NT_CONTIG {
                    if !training_seq.is_empty() {
                        training_seq.extend_from_slice(b"TTAATTAATTAA");
                    }
                    training_seq.extend_from_slice(&tmpvec);
                }
            }
            let encoded_training = if self.mask_n_runs {
                EncodedSequence::with_masking(&training_seq)
            } else {
                EncodedSequence::without_masking(&training_seq)
            };
            let mut untrained = UntrainedOrphos::with_config(analyser.config.clone())
                .expect("Failed to configure orphos for training");
            let training = untrained
                .train_single_genome(&encoded_training)
                .expect("Training failed")
                .into_training();

            log::info!("Single mode: analysis pass...");
            let mut cursor2 = Cursor::new(&self.raw_fasta_bytes);
            let mut reader2 = open_fasta(&mut cursor2);
            while let Some(record) = reader2.next() {
                let seqrec = record.expect("Invalid FASTA record");
                let (id, desc) = seqrec.id_desc().unwrap();
                let tmpdesc = desc.map(|s| s.to_string());
                let tmpid = id.to_owned();
                let tmpvec = seqrec.full_seq().to_vec();
                if tmpvec.len() < MIN_NT_CONTIG {
                    log::warn!(
                        "Contig found with less than {} nucleotides. Ignoring...",
                        MIN_NT_CONTIG
                    );
                    continue;
                }
                let tmpres = analyser
                    .analyze_sequence_bytes_with_training(&tmpvec, tmpid, tmpdesc, training.clone())
                    .expect("Error analysing FASTA record.");
                all_results.push(tmpres);
            }
        } else {
            // ── Metagenomic/anonymous mode ──
            log::info!("Creating reader and Orphos analyser...");
            let mut cursor = Cursor::new(&self.raw_fasta_bytes);
            let mut reader = open_fasta(&mut cursor);

            log::info!("Entering while loop (analysing)...");
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid FASTA record");
                let (id, desc) = seqrec.id_desc().unwrap();
                let tmpdesc = desc.map(|s| s.to_string());
                let tmpid = id.to_owned();
                let tmpvec = seqrec.full_seq().to_vec();
                if tmpvec.len() < MIN_NT_CONTIG {
                    log::warn!(
                        "Contig found with less than {} nucleotides. Ignoring...",
                        MIN_NT_CONTIG
                    );
                    continue;
                }
                let tmpres = analyser
                    .analyze_sequence_bytes(&tmpvec, tmpid, tmpdesc)
                    .expect("Error analysing FASTA record.");
                all_results.push(tmpres);
            }
        }

        log::info!("Analysis done. Saving results as attribute...");
        self.gene_count = Some(all_results.iter().map(|r| r.genes.len()).sum::<usize>());
        self.sequence_count = Some(all_results.len());
        self.results = Some(all_results);
        self.raw_fasta_bytes = Vec::new(); // free memory
    }

    pub fn get_results(&mut self, format: String) -> String {
        let output_format = match format.as_str() {
            "gbk" | "genbank" => OutputFormat::Genbank,
            "gff" => OutputFormat::Gff,
            "sco" => OutputFormat::Sco,
            "gca" => OutputFormat::Gca,
            _ => panic!("Invalid output format"),
        };

        let mut gff_bytes = Vec::new();
        for result in self.results.as_ref().expect("No results available") {
            let _ = write_results(&mut gff_bytes, result, output_format);
        }

        self.gff_bgz = bgzf_compress(&gff_bytes);
        self.gff_csi = csi_index_gff(&self.gff_bgz);

        let mut results = json::JsonValue::new_object();
        results["output_file"] = json::JsonValue::String(
            String::from_utf8(gff_bytes).expect("Unable to obtain string of the output file"),
        );
        results["gene_count"] = json::JsonValue::Number(self.gene_count.unwrap().into());
        results["sequence_count"] = json::JsonValue::Number(self.sequence_count.unwrap().into());

        results.dump()
    }

    pub fn take_fasta_bgz(&mut self) -> Vec<u8> { std::mem::take(&mut self.fasta_bgz) }
    pub fn take_fasta_fai(&mut self) -> Vec<u8> { std::mem::take(&mut self.fasta_fai) }
    pub fn take_fasta_gzi(&mut self) -> Vec<u8> { std::mem::take(&mut self.fasta_gzi) }
    pub fn take_gff_bgz(&mut self)   -> Vec<u8> { std::mem::take(&mut self.gff_bgz) }
    pub fn take_gff_csi(&mut self)   -> Vec<u8> { std::mem::take(&mut self.gff_csi) }
}

impl OrphosData {
    fn orphos_config(&self) -> OrphosConfig {
        OrphosConfig {
            metagenomic: self.metag,
            closed_ends: self.closed_ends,
            mask_n_runs: self.mask_n_runs,
            force_non_sd: self.force_non_sd,
            quiet: true,
            output_format: self.format,
            translation_table: if self.translation_table == 0 {
                None
            } else {
                Some(self.translation_table)
            },
            num_threads: None,
        }
    }
}
