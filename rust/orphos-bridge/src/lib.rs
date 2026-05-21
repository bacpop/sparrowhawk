use orphos_core::config::{OrphosConfig, OutputFormat};
use orphos_core::engine::{OrphosAnalyzer, UntrainedOrphos};
use orphos_core::output::write_results;
use orphos_core::results::OrphosResults;
use orphos_core::sequence::encoded::EncodedSequence;
use wasm_bindgen::prelude::*;
use wasm_bindgen_file_reader::WebSysFile;

mod htslib;
use htslib::{bgzf_compress, csi_index_gff, faidx_index_fasta};

pub mod fastx_wasm;
extern crate console_error_panic_hook;

use crate::fastx_wasm::open_fasta;
use flate2::read::MultiGzDecoder;
use seq_io::fasta::Record;
use std::collections::HashMap;
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
    cds_fasta: Vec<u8>,
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
            cds_fasta: Vec::new(),
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
        self.cds_fasta.clear();

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
                append_cds_fasta_records(&tmpres, &tmpvec, &mut self.cds_fasta);
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
                append_cds_fasta_records(&tmpres, &tmpvec, &mut self.cds_fasta);
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

        let gff_bytes = self.build_gff_bytes(output_format, None);
        self.store_gff_bytes(&gff_bytes);
        self.results_json(gff_bytes)
    }

    pub fn get_annotated_results(&mut self, format: String, amr_json: String) -> String {
        let output_format = match format.as_str() {
            "gbk" | "genbank" => OutputFormat::Genbank,
            "gff" => OutputFormat::Gff,
            "sco" => OutputFormat::Sco,
            "gca" => OutputFormat::Gca,
            _ => panic!("Invalid output format"),
        };
        let annotations = parse_amr_annotations(&amr_json);
        let gff_bytes = self.build_gff_bytes(output_format, Some(&annotations));
        self.store_gff_bytes(&gff_bytes);
        self.results_json(gff_bytes)
    }

    pub fn get_cds_fasta(&self) -> String {
        String::from_utf8(self.cds_fasta.clone()).expect("Unable to obtain CDS FASTA")
    }

    pub fn get_gene_metadata_json(&self) -> String {
        gene_metadata_json(self.results.as_ref().expect("No results available"))
    }

    pub fn take_fasta_bgz(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.fasta_bgz)
    }
    pub fn take_fasta_fai(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.fasta_fai)
    }
    pub fn take_fasta_gzi(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.fasta_gzi)
    }
    pub fn take_gff_bgz(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.gff_bgz)
    }
    pub fn take_gff_csi(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.gff_csi)
    }
}

impl OrphosData {
    fn build_gff_bytes(
        &self,
        output_format: OutputFormat,
        annotations: Option<&HashMap<String, AmrAnnotation>>,
    ) -> Vec<u8> {
        let mut gff_bytes = Vec::new();
        for result in self.results.as_ref().expect("No results available") {
            let _ = write_results(&mut gff_bytes, result, output_format);
        }
        if matches!(output_format, OutputFormat::Gff) {
            if let Some(annotations) = annotations {
                gff_bytes = annotate_gff_bytes(&gff_bytes, annotations);
            }
        }
        gff_bytes
    }

    fn store_gff_bytes(&mut self, gff_bytes: &[u8]) {
        self.gff_bgz = bgzf_compress(gff_bytes);
        self.gff_csi = csi_index_gff(&self.gff_bgz);
    }

    fn results_json(&self, gff_bytes: Vec<u8>) -> String {
        let mut results = json::JsonValue::new_object();
        results["output_file"] = json::JsonValue::String(
            String::from_utf8(gff_bytes).expect("Unable to obtain string of the output file"),
        );
        results["gene_count"] = json::JsonValue::Number(self.gene_count.unwrap().into());
        results["sequence_count"] = json::JsonValue::Number(self.sequence_count.unwrap().into());
        results.dump()
    }

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

#[derive(Debug, Clone)]
struct AmrAnnotation {
    unit_id: String,
    unit_label: String,
    unit_type: String,
    call_type: String,
    element_symbol: String,
    gene_symbol: String,
    allele_symbol: String,
    family: String,
    class_name: String,
    subclass: String,
    call_fraction: f64,
    diagnostic_matched: usize,
    diagnostic_total: usize,
}

impl AmrAnnotation {
    fn score(&self) -> (u8, u64, usize) {
        let specificity = if self.unit_type == "exact_gene" { 1 } else { 0 };
        (
            specificity,
            (self.call_fraction * 1_000_000.0).round() as u64,
            self.diagnostic_matched,
        )
    }
}

fn append_cds_fasta_records(result: &OrphosResults, sequence: &[u8], output: &mut Vec<u8>) {
    for gene in &result.genes {
        let begin = gene
            .coordinates
            .begin
            .min(gene.coordinates.end)
            .saturating_sub(1);
        let end = gene
            .coordinates
            .begin
            .max(gene.coordinates.end)
            .min(sequence.len());
        if begin >= end {
            continue;
        }
        output.extend_from_slice(b">");
        output.extend_from_slice(gene.annotation.identifier.as_bytes());
        output.extend_from_slice(b"\n");

        let slice = &sequence[begin..end];
        let cds = if format!("{:?}", gene.coordinates.strand) == "Reverse" {
            reverse_complement(slice)
        } else {
            slice.iter().map(|base| base.to_ascii_uppercase()).collect()
        };
        for chunk in cds.chunks(80) {
            output.extend_from_slice(chunk);
            output.extend_from_slice(b"\n");
        }
    }
}

fn gene_metadata_json(results: &[OrphosResults]) -> String {
    let mut out = json::JsonValue::new_object();
    for result in results {
        for gene in &result.genes {
            let mut metadata = json::JsonValue::new_object();
            metadata["contig"] = json::JsonValue::String(result.sequence_info.header.clone());
            metadata["start"] = json::JsonValue::Number(gene.coordinates.begin.into());
            metadata["end"] = json::JsonValue::Number(gene.coordinates.end.into());
            metadata["strand"] = json::JsonValue::String(
                match format!("{:?}", gene.coordinates.strand).as_str() {
                    "Forward" => "+",
                    "Reverse" => "-",
                    _ => ".",
                }
                .to_string(),
            );
            out[&gene.annotation.identifier] = metadata;
        }
    }
    out.dump()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            other => other,
        })
        .collect()
}

fn parse_amr_annotations(amr_json: &str) -> HashMap<String, AmrAnnotation> {
    let parsed = match json::parse(amr_json) {
        Ok(value) => value,
        Err(_) => return HashMap::new(),
    };
    let mut annotations: HashMap<String, AmrAnnotation> = HashMap::new();
    for hit in parsed["hits"].members() {
        let query_id = hit["query_id"].as_str().unwrap_or_default().to_string();
        if query_id.is_empty() {
            continue;
        }
        let annotation = AmrAnnotation {
            unit_id: json_string(hit, "unit_id"),
            unit_label: json_string(hit, "unit_label"),
            unit_type: json_string(hit, "unit_type"),
            call_type: json_string(hit, "call_type"),
            element_symbol: json_string(hit, "element_symbol"),
            gene_symbol: json_string(hit, "gene_symbol"),
            allele_symbol: json_string(hit, "allele_symbol"),
            family: json_string(hit, "family"),
            class_name: json_string(hit, "class_name"),
            subclass: json_string(hit, "subclass"),
            call_fraction: hit["call_fraction"].as_f64().unwrap_or(0.0),
            diagnostic_matched: hit["first_pass_distinct"].as_usize().unwrap_or(0),
            diagnostic_total: hit["first_pass_diagnostic_total"].as_usize().unwrap_or(0),
        };
        match annotations.get(&query_id) {
            Some(current) if current.score() >= annotation.score() => {}
            _ => {
                annotations.insert(query_id, annotation);
            }
        }
    }
    annotations
}

fn json_string(hit: &json::JsonValue, key: &str) -> String {
    hit[key].as_str().unwrap_or_default().to_string()
}

fn annotate_gff_bytes(gff_bytes: &[u8], annotations: &HashMap<String, AmrAnnotation>) -> Vec<u8> {
    let gff = String::from_utf8_lossy(gff_bytes);
    let mut output = String::with_capacity(gff.len() + annotations.len() * 256);
    for line in gff.lines() {
        if line.starts_with('#') {
            output.push_str(line);
            output.push('\n');
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() != 9 || columns[2] != "CDS" {
            output.push_str(line);
            output.push('\n');
            continue;
        }
        if let Some(gene_id) = gff_id(columns[8]) {
            if let Some(annotation) = annotations.get(&gene_id) {
                let mut attrs = columns[8].trim_end_matches(';').to_string();
                if !attrs.is_empty() {
                    attrs.push(';');
                }
                attrs.push_str(&format!(
                    "Name={};amr_unit_id={};amr_unit_label={};amr_unit_type={};amr_call_type={};amr_element_symbol={};amr_gene_symbol={};amr_allele_symbol={};amr_family={};amr_class={};amr_subclass={};amr_call_fraction={:.4};amr_diagnostic_kmers={}/{}",
                    gff_escape(&annotation.unit_label),
                    gff_escape(&annotation.unit_id),
                    gff_escape(&annotation.unit_label),
                    gff_escape(&annotation.unit_type),
                    gff_escape(&annotation.call_type),
                    gff_escape(&annotation.element_symbol),
                    gff_escape(&annotation.gene_symbol),
                    gff_escape(&annotation.allele_symbol),
                    gff_escape(&annotation.family),
                    gff_escape(&annotation.class_name),
                    gff_escape(&annotation.subclass),
                    annotation.call_fraction,
                    annotation.diagnostic_matched,
                    annotation.diagnostic_total,
                ));
                output.push_str(&columns[..8].join("\t"));
                output.push('\t');
                output.push_str(&attrs);
                output.push('\n');
                continue;
            }
        }
        output.push_str(line);
        output.push('\n');
    }
    output.into_bytes()
}

fn gff_id(attributes: &str) -> Option<String> {
    attributes
        .split(';')
        .find_map(|field| field.strip_prefix("ID=").map(|id| id.to_string()))
}

fn gff_escape(value: &str) -> String {
    let mut escaped = String::new();
    for byte in value.bytes() {
        match byte {
            b'\t' | b'\n' | b'\r' | b'%' | b';' | b'=' | b'&' | b',' | b'#' | b'?' => {
                escaped.push_str(&format!("%{byte:02X}"))
            }
            0x20..=0x7e => escaped.push(byte as char),
            other => escaped.push_str(&format!("%{other:02X}")),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_handles_lowercase_and_unknowns() {
        assert_eq!(reverse_complement(b"acgtn"), b"NACGT");
    }

    #[test]
    fn annotate_gff_adds_best_amr_hit_by_gene_id() {
        let gff = b"##gff-version  3\ncontig\tOrphos\tCDS\t1\t12\t.\t+\t0\tID=gene_1;partial=00;\n";
        let amr_json = r#"{"hits":[{"query_id":"gene_1","unit_id":"u1","unit_label":"blaA","unit_type":"hierarchy_node","call_type":"family","element_symbol":"","gene_symbol":"blaA","allele_symbol":"","family":"bla","class_name":"BETA-LACTAM","subclass":"","call_fraction":0.5,"first_pass_distinct":5,"first_pass_diagnostic_total":10,"member_count":3},{"query_id":"gene_1","unit_id":"u2","unit_label":"aac(6')-Ib","unit_type":"exact_gene","call_type":"gene","element_symbol":"aac(6')-Ib","gene_symbol":"aac","allele_symbol":"aac(6')-Ib","family":"aac","class_name":"AMINOGLYCOSIDE","subclass":"","call_fraction":0.4,"first_pass_distinct":4,"first_pass_diagnostic_total":10,"member_count":1}]}"#;
        let annotations = parse_amr_annotations(amr_json);
        let annotated = String::from_utf8(annotate_gff_bytes(gff, &annotations)).unwrap();
        assert!(annotated.contains("ID=gene_1;partial=00;Name=aac(6')-Ib;"));
        assert!(annotated.contains("amr_unit_type=exact_gene"));
        assert!(annotated.contains("amr_diagnostic_kmers=4/10"));
    }

    #[test]
    fn gff_escape_keeps_printable_gene_symbols_and_escapes_separators() {
        assert_eq!(gff_escape("aac(6')-Ib"), "aac(6')-Ib");
        assert_eq!(
            gff_escape("a;b=c&d,e%f#q?\tg"),
            "a%3Bb%3Dc%26d%2Ce%25f%23q%3F%09g"
        );
    }
}
