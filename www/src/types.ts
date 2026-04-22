export type Dict<T> = Record<string, T>

export interface WorkerState {
    worker: Worker | null;
    worker_ska: Worker | null;
    worker_deacon: Worker | null;
    worker_amr: Worker | null;
    workers_sketchlib: Worker[];
    workers_orphos: Worker[];
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

export interface AllResultsSka {
    alignResults: Dict<Alignment>
    mapResults: Dict<IsolateMapping>
    ref: string[]
    error: string | null
    clusterResults: Dict<number> | null
    transmissionGraph: TransmissionGraphData | null
}

export interface SampleIdentifyResult {
    idProbs: number[]
    idSpecies: string[]
    idMetadata: string[]
}

export interface AllResultsSketchlib {
    results: Dict<SampleIdentifyResult>
    error: string | null
}

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

export interface AmrHit {
    contig: string;
    gene_id: string | null;
    gene_family: string;
    class_name: string | null;
    start: number;
    end: number;
    distinct_hit_kmers: number;
    total_hit_kmers: number;
    reference_coverage_pct: number;
    call_type: string;
}

export interface AmrSampleResult {
    sample_name: string;
    database_profile: string;
    hits: AmrHit[];
    gene_count: number;
    family_count: number;
}

export interface AllResultsAmr {
    indexFileName: string | null;
    indexInfo: string | null;
    indexLoaded: boolean;
    isLoadingIndex: boolean;
    results: Dict<AmrSampleResult>;
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
}
