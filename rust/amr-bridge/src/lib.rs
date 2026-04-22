#[cfg(target_arch = "wasm32")]
use js_sys::{Array, Object, Reflect, Uint8Array};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneEntry {
    pub id: String,
    pub family: String,
    pub class_name: String,
    pub length: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum KmerAssignment {
    Gene(usize),
    Family(usize),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrIndex {
    pub k: usize,
    pub genes: Vec<GeneEntry>,
    pub families: Vec<String>,
    pub kmer_map: HashMap<u64, KmerAssignment>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrHit {
    pub contig: String,
    pub gene_id: Option<String>,
    pub gene_family: String,
    pub class_name: Option<String>,
    pub start: usize,
    pub end: usize,
    pub distinct_hit_kmers: usize,
    pub total_hit_kmers: usize,
    pub reference_coverage_pct: f64,
    pub call_type: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrSampleResult {
    pub sample_name: String,
    pub database_profile: String,
    pub hits: Vec<AmrHit>,
    pub gene_count: usize,
    pub family_count: usize,
}

#[derive(Debug, Clone)]
pub struct DetectParams {
    pub min_gene_hits: usize,
    pub min_family_hits: usize,
}

impl Default for DetectParams {
    fn default() -> Self {
        Self {
            min_gene_hits: 8,
            min_family_hits: 12,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct GeneAccumulator {
    count: usize,
    positions: HashSet<usize>,
    min_pos: usize,
    max_pos: usize,
}

impl GeneAccumulator {
    fn add(&mut self, pos: usize, k: usize) {
        self.count += 1;
        self.positions.insert(pos);
        if self.count == 1 {
            self.min_pos = pos;
            self.max_pos = pos + k;
        } else {
            self.min_pos = self.min_pos.min(pos);
            self.max_pos = self.max_pos.max(pos + k);
        }
    }
}

#[derive(Debug, Default, Clone)]
struct FamilyAccumulator {
    count: usize,
    positions: HashSet<usize>,
    min_pos: usize,
    max_pos: usize,
}

impl FamilyAccumulator {
    fn add(&mut self, pos: usize, k: usize) {
        self.count += 1;
        self.positions.insert(pos);
        if self.count == 1 {
            self.min_pos = pos;
            self.max_pos = pos + k;
        } else {
            self.min_pos = self.min_pos.min(pos);
            self.max_pos = self.max_pos.max(pos + k);
        }
    }
}

#[derive(Debug, Clone)]
struct FastaRecord {
    name: String,
    seq: Vec<u8>,
}

fn encode_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn canonical_kmer(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() || bytes.len() > 31 {
        return None;
    }
    let mut fwd = 0u64;
    let mut rev = 0u64;
    for &base in bytes {
        let code = encode_base(base)? as u64;
        fwd = (fwd << 2) | code;
        let comp = 3 - code;
        rev = (rev >> 2) | (comp << ((bytes.len() - 1) * 2));
    }
    Some(fwd.min(rev))
}

fn parse_fasta_bytes(bytes: &[u8]) -> Result<Vec<FastaRecord>, String> {
    let text = std::str::from_utf8(bytes).map_err(|e| format!("Invalid UTF-8 FASTA: {e}"))?;
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                records.push(FastaRecord {
                    name,
                    seq: std::mem::take(&mut current_seq),
                });
            }
            current_name = Some(
                rest.split_whitespace()
                    .next()
                    .unwrap_or("unknown")
                    .to_string(),
            );
        } else if !line.trim().is_empty() {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = current_name {
        records.push(FastaRecord {
            name,
            seq: current_seq,
        });
    }
    if records.is_empty() {
        return Err("No FASTA records found".to_string());
    }
    Ok(records)
}

fn family_from_header(header: &str) -> String {
    header.split('_').next().unwrap_or(header).to_string()
}

fn parse_config_map(root: &Path) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let cfg = root.join("config");
    let Ok(text) = fs::read_to_string(cfg) else {
        return map;
    };
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() >= 2 {
            map.insert(parts[0].to_string(), parts[1].to_string());
        }
    }
    map
}

pub fn build_index_from_resfinder_db(root: &Path, k: usize) -> Result<AmrIndex, String> {
    if k == 0 || k > 31 {
        return Err("k must be between 1 and 31".to_string());
    }
    let class_map = parse_config_map(root);
    let config_path = root.join("config");
    let config_text = fs::read_to_string(&config_path)
        .map_err(|e| format!("Failed to read {}: {e}", config_path.display()))?;

    let mut genes = Vec::new();
    let mut family_to_id = HashMap::<String, usize>::new();
    let mut families = Vec::<String>::new();
    let mut raw_kmers = HashMap::<u64, Vec<(usize, usize)>>::new();

    for line in config_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.is_empty() {
            continue;
        }
        let prefix = parts[0];
        let fasta_path = root.join(format!("{prefix}.fsa"));
        if !fasta_path.exists() {
            continue;
        }
        let class_name = class_map
            .get(prefix)
            .cloned()
            .unwrap_or_else(|| prefix.to_string());
        let bytes = fs::read(&fasta_path)
            .map_err(|e| format!("Failed to read {}: {e}", fasta_path.display()))?;
        for record in parse_fasta_bytes(&bytes)? {
            let family = family_from_header(&record.name);
            let family_id = *family_to_id.entry(family.clone()).or_insert_with(|| {
                families.push(family.clone());
                families.len() - 1
            });
            let gene_id = genes.len();
            genes.push(GeneEntry {
                id: record.name.clone(),
                family,
                class_name: class_name.clone(),
                length: record.seq.len(),
            });
            if record.seq.len() < k {
                continue;
            }
            let mut seen_gene = HashSet::new();
            for pos in 0..=record.seq.len() - k {
                if let Some(code) = canonical_kmer(&record.seq[pos..pos + k]) {
                    if seen_gene.insert(code) {
                        raw_kmers
                            .entry(code)
                            .or_default()
                            .push((gene_id, family_id));
                    }
                }
            }
        }
    }

    let mut kmer_map = HashMap::new();
    for (kmer, refs) in raw_kmers {
        let gene_ids: HashSet<usize> = refs.iter().map(|(gid, _)| *gid).collect();
        if gene_ids.len() == 1 {
            kmer_map.insert(kmer, KmerAssignment::Gene(*gene_ids.iter().next().unwrap()));
            continue;
        }
        let family_ids: HashSet<usize> = refs.iter().map(|(_, fid)| *fid).collect();
        if family_ids.len() == 1 {
            kmer_map.insert(
                kmer,
                KmerAssignment::Family(*family_ids.iter().next().unwrap()),
            );
        }
    }

    Ok(AmrIndex {
        k,
        genes,
        families,
        kmer_map,
    })
}

pub fn save_index(index: &AmrIndex, path: &Path) -> Result<(), String> {
    let bytes = bincode::serialize(index).map_err(|e| format!("Serialize index: {e}"))?;
    fs::write(path, bytes).map_err(|e| format!("Write index {}: {e}", path.display()))
}

pub fn load_index(path: &Path) -> Result<AmrIndex, String> {
    let bytes = fs::read(path).map_err(|e| format!("Read index {}: {e}", path.display()))?;
    bincode::deserialize(&bytes).map_err(|e| format!("Deserialize index: {e}"))
}

pub fn detect_fasta(
    index: &AmrIndex,
    bytes: &[u8],
    sample_name: &str,
    params: &DetectParams,
) -> Result<AmrSampleResult, String> {
    let records = parse_fasta_bytes(bytes)?;
    let mut hits = Vec::new();
    for record in records {
        let mut gene_hits = HashMap::<usize, GeneAccumulator>::new();
        let mut family_hits = HashMap::<usize, FamilyAccumulator>::new();

        if record.seq.len() >= index.k {
            for pos in 0..=record.seq.len() - index.k {
                let Some(kmer) = canonical_kmer(&record.seq[pos..pos + index.k]) else {
                    continue;
                };
                match index.kmer_map.get(&kmer) {
                    Some(KmerAssignment::Gene(gene_id)) => {
                        gene_hits.entry(*gene_id).or_default().add(pos, index.k);
                    }
                    Some(KmerAssignment::Family(family_id)) => {
                        family_hits.entry(*family_id).or_default().add(pos, index.k);
                    }
                    None => {}
                }
            }
        }

        let mut claimed_families = HashSet::new();
        for (gene_id, acc) in gene_hits {
            if acc.positions.len() < params.min_gene_hits {
                continue;
            }
            let gene = &index.genes[gene_id];
            let coverage = (acc.positions.len() as f64 * index.k as f64 / gene.length as f64
                * 100.0)
                .min(100.0);
            claimed_families.insert(gene.family.clone());
            hits.push(AmrHit {
                contig: record.name.clone(),
                gene_id: Some(gene.id.clone()),
                gene_family: gene.family.clone(),
                class_name: Some(gene.class_name.clone()),
                start: acc.min_pos,
                end: acc.max_pos,
                distinct_hit_kmers: acc.positions.len(),
                total_hit_kmers: acc.count,
                reference_coverage_pct: coverage,
                call_type: "gene".to_string(),
            });
        }

        for (family_id, acc) in family_hits {
            let family = &index.families[family_id];
            if claimed_families.contains(family) || acc.positions.len() < params.min_family_hits {
                continue;
            }
            hits.push(AmrHit {
                contig: record.name.clone(),
                gene_id: None,
                gene_family: family.clone(),
                class_name: None,
                start: acc.min_pos,
                end: acc.max_pos,
                distinct_hit_kmers: acc.positions.len(),
                total_hit_kmers: acc.count,
                reference_coverage_pct: 0.0,
                call_type: "family".to_string(),
            });
        }
    }

    let gene_count = hits.iter().filter(|hit| hit.call_type == "gene").count();
    let family_count = hits.iter().filter(|hit| hit.call_type == "family").count();
    Ok(AmrSampleResult {
        sample_name: sample_name.to_string(),
        database_profile: "resfinder".to_string(),
        hits,
        gene_count,
        family_count,
    })
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub struct WasmAmrIndex {
    inner: AmrIndex,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl WasmAmrIndex {
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<WasmAmrIndex, JsValue> {
        console_error_panic_hook::set_once();
        let inner: AmrIndex = bincode::deserialize(data)
            .map_err(|e| JsValue::from_str(&format!("Bad AMR index: {e}")))?;
        Ok(Self { inner })
    }

    pub fn info(&self) -> String {
        format!(
            "k={} ({} genes, {} families, {} kmers)",
            self.inner.k,
            self.inner.genes.len(),
            self.inner.families.len(),
            self.inner.kmer_map.len()
        )
    }
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub struct WasmAmrSession {
    index: AmrIndex,
    pending: Vec<u8>,
    min_gene_hits: usize,
    min_family_hits: usize,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl WasmAmrSession {
    #[wasm_bindgen(constructor)]
    pub fn new(
        index: &WasmAmrIndex,
        min_gene_hits: usize,
        min_family_hits: usize,
    ) -> WasmAmrSession {
        WasmAmrSession {
            index: index.inner.clone(),
            pending: Vec::new(),
            min_gene_hits,
            min_family_hits,
        }
    }

    pub fn push_chunk(&mut self, chunk: &[u8]) {
        self.pending.extend_from_slice(chunk);
    }

    pub fn finish(&mut self, sample_name: &str) -> Result<JsValue, JsValue> {
        let params = DetectParams {
            min_gene_hits: self.min_gene_hits,
            min_family_hits: self.min_family_hits,
        };
        let result = detect_fasta(&self.index, &self.pending, sample_name, &params)
            .map_err(|e| JsValue::from_str(&e))?;
        serde_to_js(&result)
    }
}

#[cfg(target_arch = "wasm32")]
fn serde_to_js<T: Serialize>(value: &T) -> Result<JsValue, JsValue> {
    let json = serde_json::to_string(value).map_err(|e| JsValue::from_str(&e.to_string()))?;
    js_sys::JSON::parse(&json).map_err(|e| e.into())
}

#[cfg(target_arch = "wasm32")]
fn js_result_summary(sample_name: &str, result: &AmrSampleResult) -> Result<JsValue, JsValue> {
    let obj = Object::new();
    set_field(&obj, "detected", true)?;
    set_field(&obj, "sampleName", sample_name)?;
    let result_js = serde_to_js(result)?;
    Reflect::set(&obj, &JsValue::from_str("result"), &result_js)?;
    Ok(obj.into())
}

#[cfg(target_arch = "wasm32")]
fn set_field<T: Into<JsValue>>(obj: &Object, key: &str, value: T) -> Result<(), JsValue> {
    Reflect::set(obj, &JsValue::from_str(key), &value.into())?;
    Ok(())
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BaselineHit {
    pub gene_id: String,
    pub family: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct EvalReport {
    pub detector_full_ids: Vec<String>,
    pub detector_hits: Vec<String>,
    pub detector_families: Vec<String>,
    pub baseline_hits: Vec<String>,
    pub baseline_families: Vec<String>,
    pub exact_overlap: Vec<String>,
    pub family_overlap: Vec<String>,
    pub detector_only: Vec<String>,
    pub baseline_only: Vec<String>,
}

pub fn parse_resfinder_json(path: &Path) -> Result<Vec<BaselineHit>, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("Read baseline {}: {e}", path.display()))?;
    let value: serde_json::Value =
        serde_json::from_str(&text).map_err(|e| format!("Parse baseline JSON: {e}"))?;
    let mut out = Vec::new();
    collect_resfinder_hits(&value, &mut out);
    out.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
    out.dedup_by(|a, b| a.gene_id == b.gene_id);
    Ok(out)
}

fn collect_resfinder_hits(value: &serde_json::Value, out: &mut Vec<BaselineHit>) {
    match value {
        serde_json::Value::Object(map) => {
            if let Some(name) = map.get("name").and_then(|v| v.as_str()) {
                let gene = map.get("gene").and_then(|v| v.as_bool()).unwrap_or(false);
                if gene {
                    out.push(BaselineHit {
                        gene_id: name.to_string(),
                        family: family_from_header(name),
                    });
                }
            }
            for child in map.values() {
                collect_resfinder_hits(child, out);
            }
        }
        serde_json::Value::Array(arr) => {
            for child in arr {
                collect_resfinder_hits(child, out);
            }
        }
        _ => {}
    }
}

pub fn compare_to_baseline(result: &AmrSampleResult, baseline: &[BaselineHit]) -> EvalReport {
    let detector_full_ids: HashSet<String> = result
        .hits
        .iter()
        .filter_map(|hit| hit.gene_id.clone())
        .collect();
    let detector_hits: HashSet<String> = detector_full_ids
        .iter()
        .map(|id| family_from_header(id))
        .collect();
    let detector_families: HashSet<String> = result
        .hits
        .iter()
        .map(|hit| hit.gene_family.clone())
        .collect();
    let baseline_hits: HashSet<String> = baseline.iter().map(|hit| hit.gene_id.clone()).collect();
    let baseline_families: HashSet<String> =
        baseline.iter().map(|hit| hit.family.clone()).collect();

    let mut exact_overlap: Vec<String> = detector_hits
        .intersection(&baseline_hits)
        .cloned()
        .collect();
    exact_overlap.sort();
    let mut family_overlap: Vec<String> = detector_families
        .intersection(&baseline_families)
        .cloned()
        .collect();
    family_overlap.sort();
    let mut detector_only: Vec<String> =
        detector_hits.difference(&baseline_hits).cloned().collect();
    detector_only.sort();
    let mut baseline_only: Vec<String> =
        baseline_hits.difference(&detector_hits).cloned().collect();
    baseline_only.sort();

    let mut detector_full_ids_vec: Vec<String> = detector_full_ids.into_iter().collect();
    detector_full_ids_vec.sort();
    let mut detector_hits_vec: Vec<String> = detector_hits.into_iter().collect();
    detector_hits_vec.sort();
    let mut detector_families_vec: Vec<String> = detector_families.into_iter().collect();
    detector_families_vec.sort();
    let mut baseline_hits_vec: Vec<String> = baseline_hits.into_iter().collect();
    baseline_hits_vec.sort();
    let mut baseline_families_vec: Vec<String> = baseline_families.into_iter().collect();
    baseline_families_vec.sort();

    EvalReport {
        detector_full_ids: detector_full_ids_vec,
        detector_hits: detector_hits_vec,
        detector_families: detector_families_vec,
        baseline_hits: baseline_hits_vec,
        baseline_families: baseline_families_vec,
        exact_overlap,
        family_overlap,
        detector_only,
        baseline_only,
    }
}

pub fn run_native_resfinder(
    resfinder_root: &Path,
    fasta: &Path,
    db_root: &Path,
    out_dir: &Path,
    blastn_path: &Path,
    kma_path: Option<&Path>,
    db_root_kma: Option<&Path>,
) -> Result<PathBuf, String> {
    fs::create_dir_all(out_dir).map_err(|e| format!("Create {}: {e}", out_dir.display()))?;
    let json_path = out_dir.join("resfinder.json");
    let module_root = resfinder_root.join("src");
    let mut cmd = std::process::Command::new("python3");
    cmd.arg("-m")
        .arg("resfinder")
        .arg("-ifa")
        .arg(fasta)
        .arg("-o")
        .arg(out_dir)
        .arg("-j")
        .arg(&json_path)
        .arg("--acquired")
        .arg("--blastPath")
        .arg(blastn_path)
        .arg("--db_path_res")
        .arg(db_root)
        .current_dir(resfinder_root)
        .env("PYTHONPATH", module_root);
    if let Some(kma_path) = kma_path {
        cmd.arg("--kmaPath").arg(kma_path);
    }
    if let Some(db_root_kma) = db_root_kma {
        cmd.arg("--db_path_res_kma").arg(db_root_kma);
    }
    let output = cmd
        .output()
        .map_err(|e| format!("Failed to run ResFinder: {e}"))?;
    if !output.status.success() {
        return Err(format!(
            "ResFinder failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    Ok(json_path)
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn detect_fasta_bytes(
    index: &WasmAmrIndex,
    data: Uint8Array,
    sample_name: &str,
) -> Result<JsValue, JsValue> {
    let mut buf = vec![0u8; data.length() as usize];
    data.copy_to(&mut buf);
    let result = detect_fasta(&index.inner, &buf, sample_name, &DetectParams::default())
        .map_err(|e| JsValue::from_str(&e))?;
    js_result_summary(sample_name, &result)
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn result_hits(result: JsValue) -> Result<Array, JsValue> {
    let hits = Reflect::get(&result, &JsValue::from_str("result"))?;
    Reflect::get(&hits, &JsValue::from_str("hits"))?.dyn_into::<Array>()
}
