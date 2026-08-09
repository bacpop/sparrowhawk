export type Dict<T> = Record<string, T>

export interface WorkerState {
    worker: Worker | null;
    worker_ska: Worker | null;
    worker_deacon: Worker | null;
    workers_sketchlib: Worker[];
    workers_orphos: Worker[];
    workers_amr: Worker[];
    workers_esm: Worker[];
}

export interface AllResults {
    nContigs: number | null;
    fastaOutput: string;
    dotOutput: string;
    gfaOutput: string;
    gfav2Output: string;
}

export interface IsolateMapping {
    mapped: boolean
    nb_variants?: number | null
    coverage?: number | null
    mapped_sequences?: string[]
    mapping_vcf?: string
}

export interface Alignment {
    aligned: boolean
    names?: string[]
    newick?: string
    alignment: string
    distances_csv?: string
}

export interface TransmissionGraphNode { id: string; cluster: number; }
export interface TransmissionGraphLink { source: string; target: string; snp_distance: number; }
export interface TransmissionGraphData {
    nodes: TransmissionGraphNode[];
    links: TransmissionGraphLink[];
}

export interface MetadataRow {
    id: string;
    date: string;
    location: string;
}

export interface TransmissionStandaloneResults {
    clusterResults: Dict<number> | null
    transmissionGraph: TransmissionGraphData | null
    error: string | null
}

export interface AllResultsSka {
    alignResults: Dict<Alignment>
    mapResults: Dict<IsolateMapping>
    mapping_vcf: string
    ref: string[]
    error: string | null
    clusterResults: Dict<number> | null
    transmissionGraph: TransmissionGraphData | null
}

export interface SampleIdentifyResult {
    idAni: number[]
    idRanks: number[]
    idSpecies: string[]
    idMetadata: string[]
}

export interface AllResultsSketchlib {
    results: Dict<SampleIdentifyResult>
    error: string | null
}

export interface GeneMetadata {
    contig: string;
    start: number;
    end: number;
    strand: string;
}

export type GeneMetadataMap = Record<string, GeneMetadata>

export interface GeneCallResult {
    fileName: string;
    outputFile: string;
    geneCount: number;
    sequenceCount: number;
    fastaBgz: Uint8Array<ArrayBuffer>;
    fastaFai: Uint8Array<ArrayBuffer>;
    fastaGzi: Uint8Array<ArrayBuffer>;
    gffBgz: Uint8Array<ArrayBuffer>;
    gffCsi: Uint8Array<ArrayBuffer>;
    geneMetadata: GeneMetadataMap;
    amrResult: AmrDetectionResult | null;
    amrHitCount: number;
    amrAnnotationTsv: string;
    amrError: string | null;
}

export interface AllResultsOrphos {
    results: Dict<GeneCallResult>;
    error: string | null
}

export interface DepletionResult {
    sampleName: string;
    totalReads: number;
    keptReads: number;
    removedReads: number;
    outputGzip: Uint8Array;
    outputGzip2: Uint8Array | null;  // null if single-end
}

export interface AllResultsDeacon {
    indexFileName: string | null;
    indexInfo: string | null;
    indexLoaded: boolean;
    isLoadingIndex: boolean;
    results: Dict<DepletionResult>;
    error: string | null
}

export interface AmrDetectionHit {
    query_id: string;
    query_kind: string;
    unit_id: string;
    unit_label: string;
    gene_id: string | null;
    element_symbol: string | null;
    gene_symbol: string | null;
    allele_symbol: string | null;
    gene_group: string;
    hierarchy_node: string | null;
    class_name: string | null;
    subclass: string | null;
    type_name: string | null;
    subtype: string | null;
    member_count: number;
    start: number;
    end: number;
    call_stage: string;
    first_pass_distinct: number;
    first_pass_total: number;
    first_pass_diagnostic_total: number;
    first_pass_fraction: number;
    call_fraction: number;
    call_type: string;
}

export interface AmrDetectionResult {
    sample_name: string;
    database_version: string;
    query_kind: string;
    index_alphabet: string;
    index_k: number;
    hits: AmrDetectionHit[];
    gene_count: number;
    gene_group_count: number;
}

export interface AllResultsAmr {
    indexFileName: string | null;
    indexInfo: string | null;
    indexLoaded: boolean;
    isLoadingIndex: boolean;
    result: AmrDetectionResult | null;
    error: string | null
}

// Minimum info needed for retry embedding, if GPU fails
export interface EsmRetry {
    file: File | null;
    sampleName: string;
}

/** One WebGPU adapter as reported by `@/platform/gpu`. */
export interface GpuAdapterInfo {
    index: number;
    name: string;
    identified: boolean; // If we managed to identify fully it (might not be the case in Firefox)
}

/** ESM-2 embeddings for one protein FASTA. */
export interface ProteinEmbeddingResult {
    sampleName: string;
    ids: string[];
    lengths: number[]; // maximum of 1022 residues, these are truncated already
    truncated: boolean[];
    dim: number;
    nSequences: number;
    backend: string;
    elapsedMs: number;
    // For setting the batches for processing with GPU
    batchMin: number;
    batchMax: number;
    // Embeddings transferred
    vectors: Float32Array;
    // Coordinates
    coords: Float32Array;
}

export interface AllResultsEsm {
    modelFileName: string | null;
    modelInfo: string | null;
    modelLoaded: boolean;
    isLoadingModel: boolean;
    modelLoadStage: string; // For reference, to see status of model loading
    backend: string | null; // Either cpu/gpu
    backendFallbackReason: string | null;
    gpuEvent: string | null; // For debugging GPU issues
    result: ProteinEmbeddingResult | null;
    error: string | null;
}

export interface ReadsPreprocessing {
    nKmers: number | null;
    histo: [];
    used_min_count: number | null;
}

export interface ProcessingState {
    isPreprocessing: boolean;
    isAssembling: boolean;
    isIndexingRef: boolean;
    isMapping: boolean;
    isMappingFiles: Set<string>;  // Track which files are being mapped
    isAligning: boolean;
    isIdentifying: boolean;
    isIdentifyingFiles: Set<string>;  // Track which files are being identified
    assemblyState: string;  // Current state from Sparrowhawk assembly
    isCallingGenes: boolean;
    isCallingGenesFiles: Set<string>;
    geneCallingStep: string;
    geneCallingProgressTotal: number;
    isFilteringDeacon: boolean;
    isFilteringDeaconFiles: Set<string>;
    isDetectingAmr: boolean;
    isDetectingAmrFiles: Set<string>;
    isClustering: boolean;
    isTransmissionStandaloneClustering: boolean;
    isEmbedding: boolean;
    isEmbeddingFiles: Set<string>;
    embeddingDone: number;
    embeddingTotal: number; 
}
