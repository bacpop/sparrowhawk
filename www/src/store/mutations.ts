import {RootState} from "@/store/state";
import type { WasmRuntimeNotice } from "@/types";
import { retainWasmRuntimeStatus } from "@/wasm/runtime";
import type { Alignment, AlignmentProgress, AlignmentStage, AlignmentStageStatus } from "@/types";
import {GeneCallResult, DepletionResult, AmrDetectionResult, Dict, TransmissionGraphData, ProteinEmbeddingResult, EsmRetry, GpuAdapterInfo} from "@/types";

export default {
    recordWasmRuntimeStatus(state: RootState, notice: WasmRuntimeNotice) {
        state.wasmRuntimeStatus[notice.moduleId] = retainWasmRuntimeStatus(
            state.wasmRuntimeStatus[notice.moduleId], notice,
        );
    },
    // Processing state mutations
    setPreprocessingState(state: RootState, isProcessing: boolean) {
        state.processingState.isPreprocessing = isProcessing;
    },
    setAssemblingState(state: RootState, isAssembling: boolean) {
        state.processingState.isAssembling = isAssembling;
    },
    setIndexingRefState(state: RootState, isIndexing: boolean) {
        state.processingState.isIndexingRef = isIndexing;
    },
    setMappingState(state: RootState, isMapping: boolean) {
        state.processingState.isMapping = isMapping;
    },
    addMappingFile(state: RootState, fileName: string) {
        state.processingState.isMappingFiles.add(fileName);
    },
    removeMappingFile(state: RootState, fileName: string) {
        console.log("Removing mapping file from list");
        state.processingState.isMappingFiles.delete(fileName);
        if (state.processingState.isMappingFiles.size === 0) {
            state.processingState.isMapping = false;
        }
    },
    setAligningState(state: RootState, isAligning: boolean) {
        state.processingState.isAligning = isAligning;
    },
    clearAlignmentLog(state: RootState) {
        state.processingState.alignmentLog = [];
    },
    recordAlignmentProgress(state: RootState, progress: AlignmentProgress) {
        const log = state.processingState.alignmentLog;
        const existingIndex = log.findIndex((entry) => entry.stage === progress.stage);
        if (existingIndex === -1) {
            log.push({ ...progress, status: "active" });
        } else {
            log.splice(existingIndex, 1, { ...progress, status: log[existingIndex].status });
        }
    },
    setAlignmentStageStatus(state: RootState, input: {
        stage: AlignmentStage;
        status: AlignmentStageStatus;
    }) {
        const entry = state.processingState.alignmentLog.find((entry) => entry.stage === input.stage);
        if (entry) entry.status = input.status;
    },
    finishAlignmentLog(state: RootState, status: "complete" | "interrupted") {
        for (const entry of state.processingState.alignmentLog) {
            if (entry.status === "active") entry.status = status;
        }
    },
    setIdentifyingState(state: RootState, isIdentifying: boolean) {
        state.processingState.isIdentifying = isIdentifying;
    },
    setAssemblyState(state: RootState, assemblyState: string) {
        state.processingState.assemblyState = assemblyState;
    },
    resetProcessingState(state: RootState) {
        state.processingState = {
            isPreprocessing: false,
            isAssembling: false,
            isIndexingRef: false,
            isMapping: false,
            isMappingFiles: new Set<string>(),
            isAligning: false,
            isObtainingAlignment: false,
            alignmentLog: [],
            isIdentifying: false,
            isIdentifyingFiles: new Set<string>(),
            assemblyState: '',
            isCallingGenes: false,
            isCallingGenesFiles: new Set<string>(),
            geneCallingStep: '',
            geneCallingProgressTotal: 0,
            isFilteringDeacon: false,
            isFilteringDeaconFiles: new Set<string>(),
            isDetectingAmr: false,
            isDetectingAmrFiles: new Set<string>(),
            isClustering: false,
            isTransmissionStandaloneClustering: false,
            isEmbedding: false,
            isEmbeddingFiles: new Set<string>(),
            embeddingDone: 0,
            embeddingTotal: 0,
        };
    },

    SET_WORKER(state: RootState, worker: Worker | null) {
        state.workerState.worker = worker;
    },
    SET_WORKER_SKA(state: RootState, worker: Worker | null) {
        state.workerState.worker_ska = worker;
    },
    SET_WORKERS_SKETCHLIB(state: RootState, workers: Worker[]) {
        state.workerState.workers_sketchlib = workers;
    },
    SET_WORKERS_ORPHOS(state: RootState, workers: Worker[]) {
        state.workerState.workers_orphos = workers;
    },
    SET_WORKER_DEACON(state: RootState, worker: Worker | null) {
        state.workerState.worker_deacon = worker;
    },
    SET_WORKERS_AMR(state: RootState, workers: Worker[]) {
        state.workerState.workers_amr = workers;
    },
    SET_WORKERS_ESM(state: RootState, workers: Worker[]) {
        state.workerState.workers_esm = workers;
    },

    setPreprocessing(state: RootState, input: { nKmers: number, histo: [], used_min_count: number, elapsedMs?: number }) {
        console.log("Preprocessing finished! Saving intermediate information in the state");
        console.log(input.nKmers);

        state.readsPreprocessing.nKmers = input.nKmers;
        state.readsPreprocessing.histo = input.histo;
        state.readsPreprocessing.used_min_count = input.used_min_count;
        state.readsPreprocessing.elapsedMs = input.elapsedMs;

    },

    setAssembly(state: RootState, input: {
        ncontigs: number,
        outfasta: string,
        outdot: string,
        outgfa: string,
        outgfav2: string,
        elapsedMs?: number,
        wasmMemoryBytes?: number
    }) {
        console.log("Assembly finished! Saving contigs as fasta in the state");

        state.allResults.nContigs = input.ncontigs;
        state.allResults.fastaOutput = input.outfasta;
        state.allResults.dotOutput = input.outdot;
        state.allResults.gfaOutput = input.outgfa;
        state.allResults.gfav2Output = input.outgfav2;
        state.allResults.elapsedMs = input.elapsedMs;
        state.allResults.wasmMemoryBytes = input.wasmMemoryBytes;

    },

    setReadsFileNames(state: RootState, input: { file1: string, file2: string }) {
        console.log("Setting names of reads");

        state.readsFileNames = input.file1 + "," + input.file2;

    },

    removeErrors(state: RootState) {
        console.log("Removing errors");

        state.errors = "";

    },

    setMemoryError(state: RootState) {
        console.log("Setting memory error in state");

        state.errors = "memory";

    },

    setFileCountError(state: RootState) {
        console.log("Setting file count error in state");

        state.errors = "file_count";

    },

    resetAllResults(state: RootState) {
        state.readsFileNames = null;
        state.min_count = 0;
        state.readsPreprocessing = {
            nKmers: null,
            histo: [],
            used_min_count: null,
        }
        state.allResults = {
            nContigs: null,
            fastaOutput: "",
            dotOutput: "",
            gfaOutput: "",
            gfav2Output: "",
        };

        if (state.workerState.worker) {
            state.workerState.worker.postMessage({reset: true});
        }
    },

    // SKA

    addRef(state: RootState, input: { name: string, sequences: string[] }) {
        console.log("vuex: Adding ref " + input.name);
        state.refSet = input.name;
        state.allResults_ska.ref = input.sequences;
    },

    addQueryFileMap(state: RootState, name: string) {
        console.log("vuex: Adding query file for mapping " + name)
        if (!state.allResults_ska.mapResults[name]) {
            state.allResults_ska.mapResults[name] = {
                mapped: true,
                nb_variants: null,
                coverage: null,
                mapped_sequences: [],
                mapping_vcf: "",
            };
        }
    },

    setMapped(state: RootState,
              input: {
                  name: string,
                  nb_variants: number | null,
                  coverage: number | null,
                  mapped_sequences: string[]
                  mapping_vcf: string
                  elapsedMs?: number
                  wasmMemoryBytes?: number
              }) {
        state.allResults_ska.mapResults[input.name].nb_variants = input.nb_variants
        state.allResults_ska.mapResults[input.name].coverage = input.coverage
        state.allResults_ska.mapResults[input.name].mapped_sequences = input.mapped_sequences
        state.allResults_ska.mapResults[input.name].mapping_vcf = input.mapping_vcf
        state.allResults_ska.mapResults[input.name].elapsedMs = input.elapsedMs
        state.allResults_ska.mapResults[input.name].wasmMemoryBytes = input.wasmMemoryBytes
        state.allResults_ska.mapping_vcf = input.mapping_vcf
    },

    setAligned(state: RootState, input: Alignment) {
        state.allResults_ska.alignResults[0] = { ...input };
        state.allResults_ska.clusterResults = null;
        state.allResults_ska.transmissionGraph = null;
    },
    invalidateAlignment(state: RootState) {
        const result = state.allResults_ska.alignResults[0];
        if (result) {
            result.aligned = false;
            result.alignmentAvailable = false;
            result.alignment_gzip = null;
            result.distances_csv_gzip = new Uint8Array();
            result.invalidated = true;
        } else {
            state.allResults_ska.alignResults[0] = {
                aligned: false, names: [], newick: "", runId: 0, k: 31, rc: false,
                alignmentAvailable: false, alignmentFrozen: false, invalidated: true,
                alignmentDownloadError: null, alignment_gzip: null, distances_csv_gzip: new Uint8Array(),
            };
        }
        state.allResults_ska.clusterResults = null;
        state.allResults_ska.transmissionGraph = null;
    },
    setObtainingAlignmentState(state: RootState, obtaining: boolean) {
        state.processingState.isObtainingAlignment = obtaining;
    },
    setAlignmentFrozen(state: RootState) {
        const result = state.allResults_ska.alignResults[0];
        if (result) result.alignmentFrozen = true;
    },
    setAlignmentDownloadError(state: RootState, message: string | null) {
        const result = state.allResults_ska.alignResults[0];
        if (result) result.alignmentDownloadError = message;
    },
    cacheAlignmentDownload(state: RootState, input: { runId: number; bytes: Uint8Array }) {
        const result = state.allResults_ska.alignResults[0];
        if (result?.runId === input.runId) result.alignment_gzip = input.bytes;
    },

    setSkaMappingError(state: RootState, msg: string) { state.allResults_ska.error = msg; },
    setSkaAlignError(state: RootState, msg: string)   { state.allResults_ska.error = msg; },

    setClusteringState(state: RootState, isClustering: boolean) {
        state.processingState.isClustering = isClustering;
    },
    setClusterResults(state: RootState, input: { clusters: Dict<number>, graph: TransmissionGraphData }) {
        state.allResults_ska.clusterResults = input.clusters;
        state.allResults_ska.transmissionGraph = input.graph;
    },
    setSkaClusterError(state: RootState, msg: string) { state.allResults_ska.error = msg; },

    setTransmissionStandaloneClusteringState(state: RootState, isClustering: boolean) {
        state.processingState.isTransmissionStandaloneClustering = isClustering;
    },
    setTransmissionStandaloneClusterResults(state: RootState, input: { clusters: Dict<number>, graph: TransmissionGraphData, distancesCapped?: boolean, elapsedMs?: number, wasmMemoryBytes?: number }) {
        state.transmissionStandalone.clusterResults = input.clusters;
        state.transmissionStandalone.transmissionGraph = input.graph;
        state.transmissionStandalone.distancesCapped = input.distancesCapped ?? false;
        state.transmissionStandalone.error = null;
        state.transmissionStandalone.elapsedMs = input.elapsedMs;
        state.transmissionStandalone.wasmMemoryBytes = input.wasmMemoryBytes;
    },
    setTransmissionStandaloneError(state: RootState, msg: string) { state.transmissionStandalone.error = msg; },
    resetTransmissionStandaloneResults(state: RootState) {
        state.transmissionStandalone = {
            clusterResults: null,
            transmissionGraph: null,
            distancesCapped: false,
            error: null,
        };
    },

    resetAllResults_ska(state: RootState) {
        state.refSet = null;
        state.allResults_ska = {
            mapResults: {},
            mapping_vcf: "",
            alignResults: {},
            ref: [],
            error: null,
            clusterResults: null,
            transmissionGraph: null,
        };
        state.processingState.isAligning = false;
        state.processingState.isObtainingAlignment = false;
        state.processingState.alignmentLog = [];

        if (state.workerState.worker_ska) {
            state.workerState.worker_ska.postMessage({reset: true});
        }
    },

    // SKETCHLIB
    addIdentifyingFile(state: RootState, sampleName: string) {
        state.processingState.isIdentifyingFiles.add(sampleName);
        state.processingState.isIdentifying = true;
    },
    removeIdentifyingFile(state: RootState, sampleName: string) {
        state.processingState.isIdentifyingFiles.delete(sampleName);
        if (state.processingState.isIdentifyingFiles.size === 0) {
            state.processingState.isIdentifying = false;
        }
    },
    saveIDResults(state: RootState, input: { sampleName: string, ani: number[], ranks: number[], names: string[], metadata: string[] }) {
        console.log("Storing results for sample: " + input.sampleName);
        state.allResults_sketchlib.results[input.sampleName] = {
            idAni: input.ani,
            idRanks: input.ranks,
            idSpecies: input.names,
            idMetadata: input.metadata,
        };
    },

    addIdentifyElapsed(state: RootState, ms: number) {
        state.allResults_sketchlib.elapsedMs = (state.allResults_sketchlib.elapsedMs ?? 0) + ms;
    },

    setIdentifyWorkerMemory(state: RootState, input: { workerIndex: number, bytes: number }) {
        state.allResults_sketchlib.memoryByWorker = {
            ...(state.allResults_sketchlib.memoryByWorker ?? {}),
            [input.workerIndex]: input.bytes,
        };
    },
    clearIdentifyWorkerMemory(state: RootState) {
        state.allResults_sketchlib.memoryByWorker = undefined;
    },

    setSketchlibError(state: RootState, msg: string) { state.allResults_sketchlib.error = msg; },

    resetAllResults_sketchlib(state: RootState) {
        state.allResults_sketchlib = {
            results: {},
            error: null,
        };

        for (const worker of state.workerState.workers_sketchlib) {
            worker.postMessage({reset: true});
        }
    },

    // ORPHOS
    addCallingGenesFile(state: RootState, fileName: string) {
        state.processingState.isCallingGenesFiles.add(fileName);
        state.processingState.isCallingGenes = true;
    },
    removeCallingGenesFile(state: RootState, fileName: string) {
        state.processingState.isCallingGenesFiles.delete(fileName);
        if (state.processingState.isCallingGenesFiles.size === 0) {
            state.processingState.isCallingGenes = false;
        }
    },
    saveGeneCallingResult(state: RootState, input: GeneCallResult) {
        state.allResults_orphos.results[input.fileName] = input;
    },

    addGeneCallingElapsed(state: RootState, ms: number) {
        state.allResults_orphos.elapsedMs = (state.allResults_orphos.elapsedMs ?? 0) + ms;
    },

    setGeneCallingWorkerMemory(state: RootState, input: { workerIndex: number, bytes: number }) {
        state.allResults_orphos.memoryByWorker = {
            ...(state.allResults_orphos.memoryByWorker ?? {}),
            [input.workerIndex]: input.bytes,
        };
    },
    clearGeneCallingWorkerMemory(state: RootState) {
        state.allResults_orphos.memoryByWorker = undefined;
    },
    setOrphosError(state: RootState, msg: string) { state.allResults_orphos.error = msg; },

    setGeneCallingStep(state: RootState, step: string) {
        state.processingState.geneCallingStep = step;
    },
    addGeneCallingProgressTotal(state: RootState, count: number) {
        state.processingState.geneCallingProgressTotal += count;
    },

    resetAllResults_orphos(state: RootState) {
        state.allResults_orphos = { results: {}, error: null };
        state.processingState.isCallingGenes = false;
        state.processingState.isCallingGenesFiles = new Set();
        state.processingState.geneCallingStep = '';
        state.processingState.geneCallingProgressTotal = 0;
        for (const worker of state.workerState.workers_orphos) {
            worker.postMessage({ reset: true });
        }
    },

    // DEACON
    setLoadingDeaconIndex(state: RootState) {
        state.allResults_deacon.isLoadingIndex = true;
    },
    setDeaconIndexLoaded(state: RootState, input: { fileName: string; info: string }) {
        state.allResults_deacon.isLoadingIndex = false;
        state.allResults_deacon.indexLoaded = true;
        state.allResults_deacon.indexFileName = input.fileName;
        state.allResults_deacon.indexInfo = input.info;
    },
    addFilteringDeaconFile(state: RootState, sampleName: string) {
        state.processingState.isFilteringDeaconFiles.add(sampleName);
        state.processingState.isFilteringDeacon = true;
    },
    removeFilteringDeaconFile(state: RootState, sampleName: string) {
        state.processingState.isFilteringDeaconFiles.delete(sampleName);
        if (state.processingState.isFilteringDeaconFiles.size === 0) {
            state.processingState.isFilteringDeacon = false;
        }
    },
    saveDeaconFilterResult(state: RootState, result: DepletionResult) {
        state.allResults_deacon.results[result.sampleName] = result;
    },
    setDeaconError(state: RootState, msg: string) { state.allResults_deacon.error = msg; },

    resetAllResults_deacon(state: RootState) {
        state.allResults_deacon = {
            indexFileName: null,
            indexInfo: null,
            indexLoaded: false,
            isLoadingIndex: false,
            results: {},
            error: null,
        };
        state.processingState.isFilteringDeacon = false;
        state.processingState.isFilteringDeaconFiles = new Set<string>();
        if (state.workerState.worker_deacon) {
            state.workerState.worker_deacon.postMessage({ reset: true });
        }
    },

    // AMR
    setLoadingAmrIndex(state: RootState) {
        state.allResults_amr.isLoadingIndex = true;
        state.allResults_amr.error = null;
    },
    setAmrIndexLoaded(state: RootState, input: { fileName: string; info: string }) {
        state.allResults_amr.isLoadingIndex = false;
        state.allResults_amr.indexLoaded = true;
        state.allResults_amr.indexFileName = input.fileName;
        state.allResults_amr.indexInfo = input.info;
    },
    addDetectingAmrFile(state: RootState, sampleName: string) {
        state.processingState.isDetectingAmrFiles.add(sampleName);
        state.processingState.isDetectingAmr = true;
        state.allResults_amr.error = null;
    },
    removeDetectingAmrFile(state: RootState, sampleName: string) {
        state.processingState.isDetectingAmrFiles.delete(sampleName);
        if (state.processingState.isDetectingAmrFiles.size === 0) {
            state.processingState.isDetectingAmr = false;
        }
    },
    saveAmrResult(state: RootState, result: AmrDetectionResult) {
        state.allResults_amr.result = result;
    },
    setAmrError(state: RootState, msg: string) {
        state.allResults_amr.isLoadingIndex = false;
        state.allResults_amr.error = msg;
    },
    resetAllResults_amr(state: RootState) {
        state.allResults_amr.result = null;
        state.allResults_amr.error = null;
        state.processingState.isDetectingAmr = false;
        state.processingState.isDetectingAmrFiles = new Set<string>();
    },

    // ESM / PROTEIN EMBEDDINGS
    setLoadingEsmModel(state: RootState, stage: string) {
        state.allResults_esm.isLoadingModel = true;
        state.allResults_esm.modelLoadStage = stage;
        state.allResults_esm.error = null;
    },
    setEsmModelLoaded(state: RootState, input: { fileName: string, info: string, backend: string }) {
        state.allResults_esm.isLoadingModel = false;
        state.allResults_esm.modelLoadStage = '';
        state.allResults_esm.modelLoaded = true;
        state.allResults_esm.modelFileName = input.fileName;
        state.allResults_esm.modelInfo = input.info;
        state.allResults_esm.backend = input.backend;
    },
   setEsmModelUnloaded(state: RootState) {
        state.allResults_esm.isLoadingModel = false;
        state.allResults_esm.modelLoadStage = '';
        state.allResults_esm.modelLoaded = false;
        state.allResults_esm.backend = null;
    },
    setEsmBackendFallback(state: RootState, reason: string) {
        state.allResults_esm.backend = "cpu";
        state.allResults_esm.backendFallbackReason = reason;
    },
    setEsmGpuEvent(state: RootState, message: string) {
        state.allResults_esm.gpuEvent = message;
    },
    setEsmRetry(state: RootState, retry: EsmRetry) {
        state.esmRetry = retry;
    },
    /** Consume the retry context, so a second abort errors instead of looping. */
    clearEsmRetry(state: RootState) {
        state.esmRetry = {file: null, sampleName: ''};
    },
    addEmbeddingFile(state: RootState, sampleName: string) {
        state.processingState.isEmbeddingFiles.add(sampleName);
        state.processingState.isEmbedding = true;
        state.allResults_esm.error = null;
    },
    removeEmbeddingFile(state: RootState, sampleName: string) {
        state.processingState.isEmbeddingFiles.delete(sampleName);
        if (state.processingState.isEmbeddingFiles.size === 0) {
            state.processingState.isEmbedding = false;
            state.processingState.embeddingDone = 0;
            state.processingState.embeddingTotal = 0;
        }
    },
    setEmbeddingProgress(state: RootState, input: { done: number, total: number }) {
        state.processingState.embeddingDone = input.done;
        state.processingState.embeddingTotal = input.total;
    },
    saveEmbeddingResult(state: RootState, result: ProteinEmbeddingResult) {
        state.allResults_esm.result = result;
    },
    setEsmError(state: RootState, msg: string) {
        state.allResults_esm.isLoadingModel = false;
        state.allResults_esm.modelLoadStage = '';
        state.allResults_esm.error = msg;
    },
    resetAllResults_esm(state: RootState) {
        state.allResults_esm.result = null;
        state.allResults_esm.error = null;
        state.processingState.isEmbedding = false;
        state.processingState.isEmbeddingFiles = new Set<string>();
        state.processingState.embeddingDone = 0;
        state.processingState.embeddingTotal = 0;
        // Note this does not remove the model loaded on the GPU, it doesn't force that for speed reasons
        for (const worker of state.workerState.workers_esm) {
            worker.postMessage({reset: true});
        }
    },

    setGpuAdapters(state: RootState, adapters: GpuAdapterInfo[]) {
        state.gpuAdapters = adapters;
    },
};
