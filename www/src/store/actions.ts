import {ActionContext} from "vuex";
import {RootState} from "@/store/state";
import WorkerSketcher from '@/workers/Sketcher.worker';
import WorkerCaller from '@/workers/Caller.worker';
import WorkerAmrDetector from '@/workers/AmrDetector.worker';
import WorkerEmbedder from '@/workers/Embedder.worker';
import WorkerAlignmentExtractor from '@/workers/AlignmentExtractor.worker';
import { observeWasmRuntime } from '@/wasm/runtime';
import {
    clampWorkerCount,
    getMaxAlignmentExtractionWorkers,
    regExpWithTwoNumbers,
    regExpForAnyFastx,
    regExpForAnyFasta,
    regExpForAnyProteinFasta,
    getFilesToProcess,
} from "@/utils";
import {listGpuAdapters as probeGpuAdapters} from "@/platform/gpu";

// GPU stall watchdog. Progress arrives per batch, so silence means the device is wedged.
const ESM_STALL_MS = 60_000;
let esmStallTimer: ReturnType<typeof setInterval> | undefined;
let esmLastProgressAt = 0;
let esmLastGpuEventAt = 0;
let esmGpuProven = false;

// Wall-clock batch timers for the pooled tabs' summary lines (main thread; batch =
// from dispatch with nothing in flight until the in-flight set empties again).
let identifyStart = 0;
let geneCallingStart = 0;
let alignmentRunCounter = 0;
let skaRequestCounter = 0;

function skaBusy(state: RootState): boolean {
    const p = state.processingState;
    return p.isAligning || p.isObtainingAlignment || p.isClustering ||
        p.isTransmissionStandaloneClustering || p.isIndexingRef || p.isMapping;
}

interface AlignmentSamplePlan {
    sampleIndex: number;
    name: string;
    firstFileIndex: number;
    secondFileIndex: number | null;
}

interface ExtractedAlignmentSample {
    type: "sampleExtracted";
    runId: number;
    sampleIndex: number;
    sampleName: string;
    packedKmers: Uint32Array;
    elapsedMs: number;
}

function clearEsmWatchdog(): void {
    clearInterval(esmStallTimer);
    esmStallTimer = undefined;
}

export default {
    async processReads(context: ActionContext<RootState, RootState>, payload: {
        acceptFiles: Array<File>,
        k: number,
        min_count: number,
        min_qual: number,
        csize: number,
        do_bloom: boolean,
        do_fit: boolean,
        no_bubble_collapse: boolean,
        no_dead_end_removal: boolean
    }) {
        const {commit, state} = context;
        console.log("Going to upload reads and assemble with k = " + payload.k + " min_count = " + payload.min_count + " min_qual = " + payload.min_qual + " csize = " + payload.csize + " do_bloom = " + payload.do_bloom)
        console.log("Checking number of uploaded files...")
        if (state.workerState.worker) {
            if (payload.acceptFiles.length > 0 && payload.acceptFiles.length < 3) {
                commit("resetAllResults");

                console.log("Removing errors if they are")
                commit("removeErrors");

                console.log("File(s) uploaded. Saving filenames...")
                commit("setReadsFileNames", {
                    file1: payload.acceptFiles[0].name,
                    file2: payload.acceptFiles[1]?.name ?? null
                });

                state.min_count = payload.min_count;
                console.log("Preprocessing...");

                // Set processing state
                commit("setPreprocessingState", true);

                state.workerState.worker.postMessage({
                    preprocess: true,
                    file1: payload.acceptFiles[0],
                    file2: payload.acceptFiles[1] ?? null,
                    min_count: payload.min_count,
                    min_qual: payload.min_qual,
                    k: payload.k,
                    verbose: true,
                    csize: payload.csize,
                    do_bloom: payload.do_bloom,
                    do_fit: payload.do_fit,
                    no_bubble_collapse: payload.no_bubble_collapse,
                    no_dead_end_removal: payload.no_dead_end_removal
                });

                state.workerState.worker.onmessage = (messageData) => {
                    // Handle state progress messages from Sparrowhawk wasm
                    if (messageData.data instanceof Object && "assemblyState" in messageData.data) {
                        console.log("[Sparrowhawk] State:", messageData.data.assemblyState);
                        commit("setAssemblyState", messageData.data.assemblyState);
                        return;
                    }

                    // Clear processing state
                    commit("setPreprocessingState", false);
                    commit("setAssemblyState", "");

                    if (messageData.data instanceof Object) {
                        if ("nKmers" in messageData.data) {
                            commit("setPreprocessing", {
                                nKmers: messageData.data.nKmers,
                                histo: messageData.data.histo,
                                used_min_count: messageData.data.used_min_count,
                                elapsedMs: messageData.data.elapsedMs,
                            });
                        } else {
                            // Something wrong has happened
                            console.log("Error found during processing, resetting all.");
                            commit("setMemoryError");
                            commit("resetAllResults");
                        }
                    }
                };
            } else {
                console.log("Wrong number of files uploaded (must be 1 or 2).");
                commit("resetAllResults");
                commit("removeErrors");
                commit("setFileCountError");
            }
        }
    },

    async doTheAssembly(context: ActionContext<RootState, RootState>) {
        const {commit, state} = context;
        console.log("Assemblying reads...")
        if (state.workerState.worker) {
            console.log("Assemblying...");

            // Set assembling state
            commit("setAssemblingState", true);

            state.workerState.worker.postMessage({
                assemble: true,
            });


            state.workerState.worker.onmessage = (messageData) => {
                // Handle state progress messages from WASM
                if ("assemblyState" in messageData.data) {
                    console.log("[Sparrowhawk] Assembly state:", messageData.data.assemblyState);
                    commit("setAssemblyState", messageData.data.assemblyState);
                    return;
                }

                // Clear assembling state
                commit("setAssemblingState", false);
                commit("setAssemblyState", "");

                commit("setAssembly", {
                    ncontigs: messageData.data.ncontigs,
                    outfasta: messageData.data.outfasta,
                    outdot: messageData.data.outdot,
                    outgfa: messageData.data.outgfa,
                    outgfav2: messageData.data.outgfav2,
                    elapsedMs: messageData.data.elapsedMs,
                    wasmMemoryBytes: messageData.data.wasmMemoryBytes,
                });
            };
        }
    },

    async resetAllResults(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
        commit("resetAllResults");
    },

    async removeErrors(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
        commit("removeErrors");
    },

    async processRef(context: ActionContext<RootState, RootState>, payload: { 
        acceptFiles: Array<File>, 
        k: number, 
        rc: boolean, 
        ambig_mask: boolean,
        repeat_mask: boolean,
        min_count: number,
        min_qual: number,
        qual_filter: number
    }) {
        const {commit, state} = context;
        if (state.processingState.isAligning || state.processingState.isObtainingAlignment ||
            state.processingState.isClustering || state.processingState.isTransmissionStandaloneClustering) return;
        console.log("Ref file uploaded, k = " + payload.k)

        // Set indexing state
        commit("setIndexingRefState", true);

        payload.acceptFiles.forEach((file: File) => {
            if (state.workerState.worker_ska) {
                state.workerState.worker_ska.postMessage({
                    ref: true, 
                    file, 
                    k: payload.k, 
                    rc: payload.rc,
                    ambig_mask: payload.ambig_mask,
                    repeat_mask: payload.repeat_mask,
                    min_count: payload.min_count,
                    min_qual: payload.min_qual,
                    qual_filter: payload.qual_filter
                });
                state.workerState.worker_ska.onmessage = (messageData) => {
                    if (!messageData.data?.ref) return;
                    // Clear indexing state
                    commit("setIndexingRefState", false);

                    console.log(messageData.data.ref.name + " has been indexed");
                    commit("addRef", {name: messageData.data.ref.name.replace(regExpForAnyFastx, ""), sequences: messageData.data.sequences});
                };
            }
        });
    },

    async processQueryMap(context: ActionContext<RootState, RootState>, payload: {
        acceptFiles: Array<File>,
        proportion_reads: number,
        min_count: number,
        min_qual: number,
        qual_filter: number
    }) {
        const {commit, state} = context;
        if (state.processingState.isAligning || state.processingState.isObtainingAlignment ||
            state.processingState.isClustering || state.processingState.isTransmissionStandaloneClustering) return;
        console.log("Query files uploaded mapping")

        // Set mapping state
        commit("setMappingState", true);
        const indxlist = getFilesToProcess(payload.acceptFiles);

        console.log(indxlist);

        indxlist.forEach((sublist: Array<number>) => {
            const messageData: any = {
                map: true,
                file: payload.acceptFiles[sublist[0]],
                revReads: (sublist.length > 1) ? payload.acceptFiles[sublist[1]] : null,
                sampleName: (sublist.length > 1) ? payload.acceptFiles[sublist[0]].name.replace(regExpWithTwoNumbers, "") : payload.acceptFiles[sublist[0]].name,
                proportion_reads: payload.proportion_reads,
                min_count: payload.min_count,
                min_qual: payload.min_qual,
                qual_filter: payload.qual_filter
            };
            // Track this file as being mapped
            commit("addMappingFile", messageData.sampleName);
            commit("addQueryFileMap", messageData.sampleName);
            if (state.workerState.worker_ska) {
                state.workerState.worker_ska.postMessage(messageData);
                state.workerState.worker_ska.onmessage = (message) => {
                    if (message.data.alignmentError || message.data.clusterError || message.data.alignmentDownloadError) return;
                    if (!("nb_variants" in message.data) && !message.data.error) return;
                    if (message.data.error) {
                        commit("setSkaMappingError", message.data.message ?? "generic");
                        commit("removeMappingFile", messageData.sampleName);
                        return;
                    }
                    console.log("Mapping variants: " + message.data.nb_variants);
                    console.log("Mapping coverage: " + message.data.coverage);

                    // Remove this file from mapping set
                    commit("removeMappingFile", message.data.name);
                    commit("setMapped", message.data);
                };
            }

        });


    },

    async processQueryAlign(context: ActionContext<RootState, RootState>, payload: {
        acceptFiles: Array<File>,
        k: number,
        proportion_reads: number,
        rc: boolean,
        min_count: number,
        min_qual: number,
        qual_filter: number,
        extraction_workers: number,
    }) {
        const {commit, state} = context;
        console.log("Processing query of uploaded files for alignment...");

        // Check first if the previous alignment has failed or it's done.
        const previous = state.allResults_ska.alignResults[0];
        if (skaBusy(state) || previous?.alignmentFrozen || previous?.invalidated) {
            commit("setSkaAlignError", "This dataset is busy or closed to additions. Clear results to start another.");
            return { success: false, invalidated: false };
        }

        // Check k values, rc
        if (previous?.aligned && (previous.k !== payload.k || previous.rc !== payload.rc)) {
            commit("setSkaAlignError", "k and reverse-complement settings must match the existing dataset.");
            return { success: false, invalidated: false };
        }
        const append = previous?.aligned === true;
        const baseIndex = append ? (previous.names?.length ?? 0) : 0;
        commit("setSkaAlignError", "");
        commit("clearAlignmentLog");

        // Get the aggregator worker
        const aggregator = state.workerState.worker_ska;
        if (!aggregator) {
            commit("setAligningState", false);
            commit("setSkaAlignError", "worker_unavailable");
            return { success: false, invalidated: false };
        }

        const runId = ++alignmentRunCounter;
        const requestedWorkers = Number.isFinite(payload.extraction_workers)
            ? Math.floor(payload.extraction_workers)
            : 2;
        const workerCount = clampWorkerCount(
            requestedWorkers,
            getMaxAlignmentExtractionWorkers(),
        );

        // Create pool (array) of workers, catching errors always if possible
        const extractors: Worker[] = [];
        try {
            const poolSize = Math.min(workerCount, Math.max(1, payload.acceptFiles.length));
            for (let index = 0; index < poolSize; index++) {
                const worker = new WorkerAlignmentExtractor();
                extractors.push(worker);
                observeWasmRuntime(worker, (notice) => commit("recordWasmRuntimeStatus", notice));
            }
        } catch (error) {
            extractors.forEach((worker) => worker.terminate());
            const message = error instanceof Error ? error.message : String(error);
            commit("setAligningState", false);
            commit("setSkaAlignError", message);
            console.error("[ska-align] Could not start extraction workers", { runId, error });
            return { success: false, invalidated: false };
        }



        // This big promise is the alignment!!
        commit("setAligningState", true);
        const pool = extractors;
        return new Promise<{ success: boolean; invalidated: boolean }>((resolve) => {
            // Inits...
            let finished = false;
            let finishSent = false;
            let samples: AlignmentSamplePlan[] | null = null;
            let nextSampleToSchedule = 0;
            let nextSampleToAdd = 0;
            let extractedCount = 0;
            let acknowledgedCount = 0;
            let poolTerminated = false;

            const runningByWorker = new Map<number, number>();
            const pendingResults = new Map<number, ExtractedAlignmentSample>();
            const acknowledgedSamples = new Set<number>();
            const wordWidth = payload.k < 32 ? 2 : 4;

            // Helper methods/functions. TODO: move out  from here (although we'd need to fix the arguments)

            const terminatePool = () => {
                if (poolTerminated) return;
                poolTerminated = true;
                pool.forEach((extractor) => extractor.terminate());
            };

            // Proper end of the promise
            const finish = (errorMessage?: string, result?: Record<string, unknown>, validationOnly = false) => {
                if (finished) return;
                finished = true;
                terminatePool();
                aggregator.removeEventListener("message", onAggregatorMessage);
                aggregator.removeEventListener("error", onAggregatorError);

                if (errorMessage) {
                    if (!validationOnly) {
                        try {
                            aggregator.postMessage({ cancelAlignment: true, runId });
                        } catch (error) {
                            console.error("[ska-align] Could not notify the worker of cancellation", error);
                        }
                        commit("invalidateAlignment");
                    }
                    console.error("[ska-align] Alignment failed", { runId, error: errorMessage });
                    commit("setSkaAlignError", errorMessage);
                } else if (result) {
                    commit("setAligned", result);
                }
                commit("finishAlignmentLog", errorMessage ? "interrupted" : "complete");
                commit("setAligningState", false);
                resolve({ success: !errorMessage, invalidated: !!errorMessage && !validationOnly });
            };

            // Properly ending, but with an error
            const fail = (error: unknown) => {
                const message = error instanceof Error ? error.message : String(error);
                finish(message || "Alignment failed.");
            };

            const maybeFinishAlignment = () => {
                if (finished || finishSent || samples === null) return;
                if (acknowledgedCount !== samples.length) return;
                finishSent = true;
                if (samples.length === 0) terminatePool();
                try {
                    aggregator.postMessage({ finishAlignment: true, runId });
                } catch (error) {
                    fail(error);
                }
            };

            // This is the method that moves the results from one worker that processed it, to the aggregator
            const flushResults = () => {
                if (finished || samples === null) return;
                while (pendingResults.has(baseIndex + nextSampleToAdd)) {
                    const sample = pendingResults.get(baseIndex + nextSampleToAdd)!;
                    pendingResults.delete(baseIndex + nextSampleToAdd);
                    try {
                        aggregator.postMessage({
                            addAlignmentSample: true,
                            runId,
                            sampleIndex: sample.sampleIndex,
                            sampleName: sample.sampleName,
                            packedKmers: sample.packedKmers,
                        }, [sample.packedKmers.buffer as ArrayBuffer]);
                    } catch (error) {
                        fail(error);
                        return;
                    }
                    nextSampleToAdd += 1;
                }
                maybeFinishAlignment();
            };

            // This manages the buffers we have here (essentially "samples" and "pendingResults"), alongside flushResults.
            const pump = () => {
                if (finished || samples === null) return;
                while (
                    nextSampleToSchedule < samples.length &&
                    nextSampleToSchedule - acknowledgedCount < pool.length
                ) {
                    const workerIndex = pool.findIndex((_, index) => !runningByWorker.has(index));
                    if (workerIndex < 0) break;

                    const sample = samples[nextSampleToSchedule];
                    const firstFile = payload.acceptFiles[sample.firstFileIndex];
                    const secondFile = sample.secondFileIndex === null
                        ? null
                        : payload.acceptFiles[sample.secondFileIndex];
                    if (!firstFile || (sample.secondFileIndex !== null && !secondFile)) {
                        fail(new Error(`The sample plan refers to a missing file for ${sample.name}.`));
                        return;
                    }

                    nextSampleToSchedule += 1;
                    runningByWorker.set(workerIndex, sample.sampleIndex);
                    try {
                        pool[workerIndex].postMessage({
                            extractSample: true,
                            runId,
                            sampleIndex: sample.sampleIndex,
                            sampleName: sample.name,
                            firstFile,
                            secondFile,
                            k: payload.k,
                            rc: payload.rc,
                            proportion_reads: payload.proportion_reads,
                            min_count: payload.min_count,
                            min_qual: payload.min_qual,
                            qual_filter: payload.qual_filter,
                        });
                    } catch (error) {
                        fail(error);
                        return;
                    }
                }
                flushResults();
            };

            const onAggregatorMessage = (event: MessageEvent) => {
                const data = event.data;
                if (!(data instanceof Object)) return;
                if (data.alignmentProgress) {
                    if (data.alignmentProgress.runId !== runId) return;
                    if (data.alignmentProgress.stage === "loading" || data.alignmentProgress.stage === "preparing") return;
                    commit("recordAlignmentProgress", data.alignmentProgress);
                    const progress = data.alignmentProgress;
                    // Finalisation stages run sequentially in the alignment worker.
                    if (progress.stage === "tree") {
                        commit("setAlignmentStageStatus", { stage: "tree-distances", status: "complete" });
                    } else if (progress.stage === "exporting") {
                        commit("setAlignmentStageStatus", { stage: "tree", status: "complete" });
                    }
                    return;
                }
                if (data.runId !== runId) return;


                // After organising all the files etc.
                if (data.alignmentPlan) {
                    if (!Array.isArray(data.samples)) {
                        fail(new Error("The alignment worker returned an invalid sample plan."));
                        return;
                    }
                    samples = data.samples as AlignmentSamplePlan[];
                    commit("recordAlignmentProgress", {
                        stage: "extracting",
                        sampleIndex: baseIndex,
                        sampleTotal: baseIndex + samples.length,
                        sampleName: null,
                        runId,
                    });
                    if (samples.length === 0) {
                        commit("setAlignmentStageStatus", { stage: "extracting", status: "complete" });
                        terminatePool();
                    }
                    pump();
                    return;
                }


                if (data.alignmentSampleAdded) {
                    if (!acknowledgedSamples.has(data.sampleIndex)) {
                        acknowledgedSamples.add(data.sampleIndex);
                        acknowledgedCount += 1;
                    }
                    if (samples !== null && acknowledgedCount === samples.length) {
                        commit("setAlignmentStageStatus", { stage: "distances", status: "complete" });
                    }
                    pump();
                    return;
                }


                if (data.aligned === true) {
                    if (!(data.distances_csv_gzip instanceof Uint8Array)) {
                        fail(new Error("The alignment worker returned invalid compressed outputs."));
                        return;
                    }
                    finish(undefined, { ...data, invalidated: false, alignmentDownloadError: null });
                    return;
                }

                if (data.error && data.alignmentError) {
                    finish(String(data.message ?? "Alignment worker failed."), undefined, data.validationOnly === true);
                }
            };

            const onAggregatorError = (event: ErrorEvent) => {
                fail(new Error(event.message || "Alignment worker failed."));
            };

            // This just defines what to do for all the workers. This should be simplified
            pool.forEach((extractor, workerIndex) => {
                extractor.onmessage = (event: MessageEvent) => {
                    const data = event.data;
                    if (!(data instanceof Object) || data.runId !== runId) return;
                    const expectedSample = runningByWorker.get(workerIndex);
                    if (data.type === "sampleError") {
                        fail(new Error(`Sample ${data.sampleName ?? expectedSample} failed: ${data.message ?? "extraction error"}`));
                        return;
                    }
                    if (data.type !== "sampleExtracted") return;
                    if (expectedSample !== data.sampleIndex) {
                        fail(new Error("An extraction worker returned a result for the wrong sample."));
                        return;
                    }

                    runningByWorker.delete(workerIndex);
                    const extracted = data as ExtractedAlignmentSample;
                    if (extracted.packedKmers.length % wordWidth !== 0) {
                        fail(new Error(`Sample ${extracted.sampleName} returned malformed packed k-mers.`));
                        return;
                    }
                    pendingResults.set(extracted.sampleIndex, extracted);
                    extractedCount += 1;
                    console.info("[ska-align] Extracted sample", {
                        runId,
                        sample: extracted.sampleName,
                        sampleIndex: extracted.sampleIndex + 1,
                        sampleTotal: baseIndex + (samples?.length ?? 0),
                        kmerCount: extracted.packedKmers.length / wordWidth,
                        elapsedMs: extracted.elapsedMs,
                    });
                    commit("recordAlignmentProgress", {
                        stage: "extracting",
                        sampleIndex: baseIndex + extractedCount,
                        sampleTotal: baseIndex + (samples?.length ?? 0),
                        sampleName: extracted.sampleName,
                        runId,
                    });
                    if (samples !== null && extractedCount === samples.length) {
                        commit("setAlignmentStageStatus", { stage: "extracting", status: "complete" });
                        terminatePool();
                    }
                    flushResults();
                    pump();
                };
                extractor.onerror = (event: ErrorEvent) => {
                    fail(new Error(`Extraction worker ${workerIndex + 1} failed: ${event.message || "worker error"}`));
                };
            });

            aggregator.addEventListener("message", onAggregatorMessage);
            aggregator.addEventListener("error", onAggregatorError);


            // And finally, the actual execution, which is done from the aggregator worker.
            try {
                aggregator.postMessage({
                    alignStart: true,
                    runId,
                    fileNames: payload.acceptFiles.map((file) => file.name),
                    k: payload.k,
                    rc: payload.rc,
                    append,
                });
            } catch (error) {
                fail(error);
                return;
            }
            console.info("[ska-align] Started extraction pool", {
                runId,
                workers: pool.length,
                sampleFiles: payload.acceptFiles.length,
                k: payload.k,
            });
        }); // Big promise
    },

    async getAlignmentDownload(context: ActionContext<RootState, RootState>): Promise<Uint8Array> {
        const { state, commit } = context;
        const result = state.allResults_ska.alignResults[0];
        const worker = state.workerState.worker_ska;

        if (skaBusy(state)) throw new Error("Wait for the current SKA operation to finish.");
        if (!worker || !result?.aligned || !result.alignmentAvailable) throw new Error("No alignment is available for download.");

        console.info("[ska-align] Obtaining alignment...");

        commit("setAlignmentFrozen");
        commit("setAlignmentDownloadError", null);
        commit("recordAlignmentProgress", { stage: "obtaining-alignment", sampleIndex: null, sampleTotal: null, sampleName: null, runId: result.runId });
        commit("setAlignmentStageStatus", { stage: "obtaining-alignment", status: "active" });

        // In case we've got it already...
        if (result.alignment_gzip) {
            commit("setAlignmentStageStatus", { stage: "obtaining-alignment", status: "complete" });
            return result.alignment_gzip;
        }
        const runId = result.runId;
        const requestId = ++skaRequestCounter;

        commit("setObtainingAlignmentState", true);
        return new Promise<Uint8Array>((resolve, reject) => {

            const finish = (bytes?: Uint8Array, error?: string) => {
                worker.removeEventListener("message", onMessage);
                worker.removeEventListener("error", onError);
                commit("setObtainingAlignmentState", false);
                commit("setAlignmentStageStatus", { stage: "obtaining-alignment", status: error ? "interrupted" : "complete" });

                if (error) { commit("setAlignmentDownloadError", error); reject(new Error(error)); }
                else { commit("cacheAlignmentDownload", { runId, bytes }); resolve(bytes!); }
            };

            const onMessage = (event: MessageEvent) => {
                const data = event.data;
                if (data?.alignmentProgress?.runId === runId && data.alignmentProgress.stage === "obtaining-alignment") {
                    commit("recordAlignmentProgress", data.alignmentProgress);
                    return;
                }
                if (data?.runId !== runId || data.requestId !== requestId) return;
                if (data.alignmentDownloadError) finish(undefined, String(data.message));
                else if (data.alignmentDownload) {
                    if (!(data.bytes instanceof Uint8Array)) finish(undefined, "Invalid alignment download buffer.");
                    else finish(data.bytes);
                }
            };

            const onError = (event: ErrorEvent) => finish(undefined, event.message || "Alignment export worker failed.");

            worker.addEventListener("message", onMessage);
            worker.addEventListener("error", onError);
            try { worker.postMessage({ exportAlignment: true, runId, requestId }); }
            catch (error) { finish(undefined, error instanceof Error ? error.message : String(error)); }
        });
    },

    async processCluster(context: ActionContext<RootState, RootState>, payload: { snp_threshold: number }) {
        const { commit, state } = context;
        if (skaBusy(state)) return;
        const worker = state.workerState.worker_ska;
        const result = state.allResults_ska.alignResults[0];
        if (!worker || !result?.aligned) { commit("setSkaClusterError", "No completed alignment dataset is available."); return; }
        const requestId = ++skaRequestCounter;
        commit("setClusteringState", true);
        return new Promise<void>((resolve) => {
            const finish = (error?: string, data?: { clusters: unknown; graph: unknown }) => {
                worker.removeEventListener("message", onMessage);
                worker.removeEventListener("error", onError);
                commit("setClusteringState", false);
                if (error) commit("setSkaClusterError", error);
                else if (data) commit("setClusterResults", data);
                resolve();
            };
            const onMessage = (event: MessageEvent) => {
                const data = event.data;
                if (data?.requestId !== requestId || data.runId !== result.runId) return;
                if (data.clusterError) finish(String(data.message));
                else if (data.clustered) finish(undefined, data);
            };
            const onError = (event: ErrorEvent) => finish(event.message || "Clustering worker failed.");
            worker.addEventListener("message", onMessage);
            worker.addEventListener("error", onError);
            try { worker.postMessage({ cluster: true, snp_threshold: payload.snp_threshold, runId: result.runId, requestId }); }
            catch (error) { finish(error instanceof Error ? error.message : String(error)); }
        });
    },

    async processTransmissionStandaloneCluster(context: ActionContext<RootState, RootState>, payload: { file: File, snp_threshold: number }) {
        const { commit, state } = context;
        if (skaBusy(state)) { commit("setTransmissionStandaloneError", "Wait for the current SKA operation to finish."); return; }
        const worker = state.workerState.worker_ska;
        if (!worker) { commit("setTransmissionStandaloneError", "worker_unavailable"); return; }
        const requestId = ++skaRequestCounter;
        commit("setTransmissionStandaloneClusteringState", true);
        commit("resetTransmissionStandaloneResults");
        return new Promise<void>((resolve) => {
            const finish = (error?: string, data?: Record<string, unknown>) => {
                worker.removeEventListener("message", onMessage);
                worker.removeEventListener("error", onError);
                commit("setTransmissionStandaloneClusteringState", false);
                if (error) commit("setTransmissionStandaloneError", error);
                else if (data) commit("setTransmissionStandaloneClusterResults", data);
                resolve();
            };
            const onMessage = (event: MessageEvent) => {
                const data = event.data;
                if (data?.requestId !== requestId || !data.standalone) return;
                if (data.clusterError) finish(String(data.message));
                else if (data.clustered) finish(undefined, data);
            };
            const onError = (event: ErrorEvent) => finish(event.message || "Clustering worker failed.");
            worker.addEventListener("message", onMessage);
            worker.addEventListener("error", onError);
            try { worker.postMessage({ transmission_cluster: true, file: payload.file, snp_threshold: payload.snp_threshold, requestId }); }
            catch (error) { finish(error instanceof Error ? error.message : String(error)); }
        });
    },

    async resetTransmissionStandaloneResults(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
        commit("resetTransmissionStandaloneResults");
    },

    async resetAllResults_ska(context: ActionContext<RootState, RootState>) {
        const {commit, state} = context;
        if (skaBusy(state)) return;
        commit("resetAllResults_ska");
    },

    async identifyFiles(context: ActionContext<RootState, RootState>, payload: { acceptFiles: Array<File>, proportion_reads: number, min_count: number, min_qual: number }) {
        const {commit, state} = context;
        console.log("Uploaded file(s) for taxonomic identification");

        const pool = state.workerState.workers_sketchlib;
        if (!pool.length) {
            console.log("No sketchlib workers available");
            return;
        }

        const indxlist = getFilesToProcess(payload.acceptFiles);

        if (indxlist.length === 0) {
            console.log("No valid samples found to process");
            return;
        }

        // A fresh batch (nothing in flight) restarts the wall clock; appended files keep it.
        if (state.processingState.isIdentifyingFiles.size === 0) {
            identifyStart = performance.now();
        }

        // Attach a result handler to each worker in the pool
        pool.forEach((worker, workerIndex) => {
            worker.onmessage = (messageData) => {
                if (messageData.data instanceof Object) {
                    if (messageData.data.error) {
                        const sampleName = messageData.data.sampleName;
                        commit("setSketchlibError", messageData.data.message ?? "generic");
                        if (sampleName) {
                            commit("removeIdentifyingFile", sampleName);
                            if (state.processingState.isIdentifyingFiles.size === 0) {
                                commit("addIdentifyElapsed", Math.round(performance.now() - identifyStart));
                            }
                        }
                    } else if ("ani" in messageData.data && "sampleName" in messageData.data) {
                        const sampleName = messageData.data.sampleName;
                        console.log("Saving results for sample: " + sampleName);
                        commit("removeIdentifyingFile", sampleName);
                        commit("saveIDResults", {
                            sampleName,
                            ani: messageData.data.ani,
                            ranks: messageData.data.ranks,
                            names: messageData.data.names,
                            metadata: messageData.data.metadata,
                        });
                        if (typeof messageData.data.wasmMemoryBytes === "number") {
                            commit("setIdentifyWorkerMemory", { workerIndex, bytes: messageData.data.wasmMemoryBytes });
                        }
                        if (state.processingState.isIdentifyingFiles.size === 0) {
                            commit("addIdentifyElapsed", Math.round(performance.now() - identifyStart));
                        }
                    } else {
                        console.log("Error found during processing");
                    }
                }
            };
        });

        // Distribute samples across the pool via round-robin
        indxlist.forEach((sublist: number[], index: number) => {
            const file1 = payload.acceptFiles[sublist[0]];
            const file2 = sublist.length > 1 ? payload.acceptFiles[sublist[1]] : null;
            const sampleName = sublist.length > 1
                ? file1.name.replace(regExpWithTwoNumbers, "")
                : file1.name.replace(regExpForAnyFastx, '');
            const worker = pool[index % pool.length];
            console.log(`Queuing identification for sample: ${sampleName} on worker ${index % pool.length}`);
            commit("addIdentifyingFile", sampleName);
            worker.postMessage({
                identify: true,
                file1,
                file2,
                sampleName,
                proportion_reads: payload.proportion_reads,
                min_count: payload.min_count,
                min_qual: payload.min_qual,
            });
        });
    },

    async initSketchlibWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
        const {commit, state} = context;
        const workerCount = clampWorkerCount(numWorkers);
        // Terminate existing workers
        for (const worker of state.workerState.workers_sketchlib) {
            worker.terminate();
        }
        commit("clearIdentifyWorkerMemory");
        // Spawn new pool
        const pool: Worker[] = [];
        for (let i = 0; i < workerCount; i++) {
            pool.push(new WorkerSketcher());
        }
        console.log(`Spawned ${workerCount} sketchlib worker(s)`);
        commit("SET_WORKERS_SKETCHLIB", pool);
    },

    async resetAllResults_sketchlib(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
        commit("resetAllResults_sketchlib");
    },

    // ORPHOS
    async initCallerWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
        const { commit, state } = context;
        const workerCount = clampWorkerCount(numWorkers);
        for (const worker of state.workerState.workers_orphos) {
            worker.terminate();
        }
        commit("clearGeneCallingWorkerMemory");
        const pool: Worker[] = [];
        for (let i = 0; i < workerCount; i++) {
            pool.push(new WorkerCaller());
        }
        pool.forEach((worker, workerIndex) => {
            worker.onmessage = (msg) => {
                if (msg.data?.error) {
                    commit("setOrphosError", msg.data.message ?? "generic");
                    if (msg.data.fileName) {
                        commit("removeCallingGenesFile", msg.data.fileName);
                        if (state.processingState.isCallingGenesFiles.size === 0) {
                            commit("addGeneCallingElapsed", Math.round(performance.now() - geneCallingStart));
                        }
                    }
                } else if (msg.data?.geneCallingStep !== undefined) {
                    commit("setGeneCallingStep", msg.data.geneCallingStep);
                } else if (msg.data?.output_file !== undefined) {
                    commit("removeCallingGenesFile", msg.data.fileName);
                    commit("saveGeneCallingResult", {
                        fileName:      msg.data.fileName,
                        outputFile:    msg.data.output_file,
                        geneCount:     msg.data.gene_count,
                        sequenceCount: msg.data.sequence_count,
                        fastaBgz:      msg.data.fasta_bgz,
                        fastaFai:      msg.data.fasta_fai,
                        fastaGzi:      msg.data.fasta_gzi,
                        gffBgz:        msg.data.gff_bgz,
                        gffCsi:        msg.data.gff_csi,
                        geneMetadata:  msg.data.gene_metadata ?? {},
                        amrResult:     msg.data.amr_result,
                        amrHitCount:   msg.data.amr_hit_count,
                        amrAnnotationTsv: msg.data.amr_tsv,
                        amrError:      msg.data.amr_error,
                    });
                    if (typeof msg.data.wasm_memory_bytes === "number") {
                        commit("setGeneCallingWorkerMemory", { workerIndex, bytes: msg.data.wasm_memory_bytes });
                    }
                    if (state.processingState.isCallingGenesFiles.size === 0) {
                        commit("addGeneCallingElapsed", Math.round(performance.now() - geneCallingStart));
                    }
                }
            };
        });
        console.log(`Spawned ${workerCount} caller worker(s)`);
        commit("SET_WORKERS_ORPHOS", pool);
    },

    async callGenes(context: ActionContext<RootState, RootState>,
        payload: {
            acceptFiles: Array<File>,
            metag: boolean,
            closed_ends: boolean,
            mask: boolean,
            tt: number,
            non_sd: boolean
            min_gene_fraction: number
            min_gene_group_fraction: number
        }) {
        const { commit, state } = context;
        const pool = state.workerState.workers_orphos;
        if (!pool.length) {
            console.log("No caller workers available");
            return;
        }
        const indxlist = getFilesToProcess(payload.acceptFiles);
        // A fresh batch (nothing in flight) restarts the wall clock; appended files keep it.
        if (state.processingState.isCallingGenesFiles.size === 0) {
            geneCallingStart = performance.now();
        }
        commit("addGeneCallingProgressTotal", indxlist.length);
        indxlist.forEach((sublist: number[], index: number) => {
            const file = payload.acceptFiles[sublist[0]];
            const fileName = file.name;
            console.log(`Queuing gene calling for: ${fileName} on worker ${index % pool.length}`);
            commit("addCallingGenesFile", fileName);
            pool[index % pool.length].postMessage({
                call: true,
                fileName,
                input_file: file,
                metag: payload.metag,
                closed_ends: payload.closed_ends,
                mask: payload.mask,
                tt: payload.tt,
                non_sd: payload.non_sd,
                min_gene_fraction: payload.min_gene_fraction,
                min_gene_group_fraction: payload.min_gene_group_fraction,
            });
        });
    },

    async resetAllResults_orphos(context: ActionContext<RootState, RootState>) {
        const { commit } = context;
        commit("resetAllResults_orphos");
    },

    // DEACON
    async loadDeaconIndex(context: ActionContext<RootState, RootState>, payload: { file: File }) {
        const { commit, state } = context;
        if (!state.workerState.worker_deacon) return;
        commit("setLoadingDeaconIndex");
        state.workerState.worker_deacon.postMessage({ loadIndex: true, file: payload.file });
        state.workerState.worker_deacon.onmessage = (msg) => {
            if (msg.data.error) {
                state.allResults_deacon.isLoadingIndex = false;
                commit("setDeaconError", msg.data.message ?? "index");
            } else if (msg.data.indexLoaded) {
                commit("setDeaconIndexLoaded", { fileName: msg.data.fileName, info: msg.data.info });
            }
        };
    },

    async filterDeaconReads(context: ActionContext<RootState, RootState>, payload: { files: Array<File>; deplete: boolean; abs_threshold: number; rel_threshold: number }) {
        const { commit, state } = context;
        if (!state.workerState.worker_deacon) return;
        const workerDeacon = state.workerState.worker_deacon;
        const { files, deplete, abs_threshold, rel_threshold } = payload;

        // Set handler once (safe to overwrite — idempotent)
        workerDeacon.onmessage = (msg) => {
            if (msg.data.error) {
                commit("setDeaconError", msg.data.message ?? "generic");
                if (msg.data.sampleName) commit("removeFilteringDeaconFile", msg.data.sampleName);
            } else if (msg.data.filtered) {
                commit("removeFilteringDeaconFile", msg.data.sampleName);
                commit("saveDeaconFilterResult", {
                    sampleName:   msg.data.sampleName,
                    totalReads:   msg.data.total,
                    keptReads:    msg.data.kept,
                    removedReads: msg.data.removed,
                    outputGzip:   msg.data.outputGzip,
                    outputGzip2:  msg.data.outputGzip2 ?? null,
                    elapsedMs:    msg.data.elapsedMs,
                    wasmMemoryBytes: msg.data.wasmMemoryBytes,
                });
            }
        };

        const indxlist = getFilesToProcess(files);
        for (const sublist of indxlist) {
            const file = files[sublist[0]];
            const revReads = sublist.length > 1 ? files[sublist[1]] : null;
            const sampleName = sublist.length > 1
                ? file.name.replace(regExpWithTwoNumbers, '')
                : file.name.replace(regExpForAnyFastx, '');

            commit("addFilteringDeaconFile", sampleName);
            workerDeacon.postMessage({ filter: true, file, revReads, deplete, abs_threshold, rel_threshold, sampleName });
        }
    },

    async resetAllResults_deacon(context: ActionContext<RootState, RootState>) {
        context.commit("resetAllResults_deacon");
    },

    // AMR
    async initAmrWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
        const {commit, state} = context;
        const workerCount = clampWorkerCount(numWorkers);
        for (const worker of state.workerState.workers_amr) {
            worker.terminate();
        }
        const pool: Worker[] = [];
        for (let i = 0; i < workerCount; i++) {
            pool.push(new WorkerAmrDetector());
        }
        console.log(`Spawned ${workerCount} AMR worker(s)`);
        commit("SET_WORKERS_AMR", pool);
    },

    async detectAmrFile(context: ActionContext<RootState, RootState>, payload: { files: Array<File>; min_gene_fraction: number; min_gene_group_fraction: number }) {
        const {commit, state} = context;
        commit("resetAllResults_amr");
        const pool = state.workerState.workers_amr;
        if (!pool.length) {
            console.log("No AMR worker available");
            commit("setAmrError", "worker");
            return;
        }
        if (payload.files.length !== 1) {
            commit("setAmrError", "file_count");
            return;
        }

        const file = payload.files[0];
        if (!regExpForAnyFasta.test(file.name)) {
            commit("setAmrError", "fasta");
            return;
        }

        const sampleName = file.name.replace(regExpForAnyFasta, "");
        const worker = pool[0];
        worker.onmessage = (messageData) => {
            if (!(messageData.data instanceof Object)) return;
            if (messageData.data.error) {
                if (messageData.data.sampleName) {
                    commit("removeDetectingAmrFile", messageData.data.sampleName);
                }
                commit("setAmrError", messageData.data.message ?? "generic");
            } else if (messageData.data.indexLoaded) {
                commit("setAmrIndexLoaded", {
                    fileName: messageData.data.fileName,
                    info: messageData.data.info,
                });
            } else if (messageData.data.detected) {
                commit("removeDetectingAmrFile", messageData.data.sampleName);
                commit("saveAmrResult", messageData.data.result);
            }
        };

        if (!state.allResults_amr.indexLoaded && !state.allResults_amr.isLoadingIndex) {
            commit("setLoadingAmrIndex");
        }
        commit("addDetectingAmrFile", sampleName);
        worker.postMessage({
            detectAmr: true,
            file,
            sampleName,
            min_gene_fraction: payload.min_gene_fraction,
            min_gene_group_fraction: payload.min_gene_group_fraction,
        });
    },

    async resetAllResults_amr(context: ActionContext<RootState, RootState>) {
        context.commit("resetAllResults_amr");
    },

    // ESM / PROTEIN EMBEDDINGS
    async initEmbedderWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
        const {commit, state, dispatch} = context;
        const workerCount = clampWorkerCount(numWorkers);
        for (const worker of state.workerState.workers_esm) {
            worker.terminate();
        }
        // The new instance holds no model, and cubecl's one-shot GPU init is reset with it.
        commit("setEsmModelUnloaded");
        const pool: Worker[] = [];
        for (let i = 0; i < workerCount; i++) {
            pool.push(new WorkerEmbedder());
        }
        
        pool.forEach((worker) => {
            worker.onmessage = (messageData) => {
                if (!(messageData.data instanceof Object)) return;
                const d = messageData.data;
                if (d.error) {
                    clearEsmWatchdog();
                    if (d.sampleName) commit("removeEmbeddingFile", d.sampleName);
                    if (d.wasmPanic) {
                        // The wasm instance has aborted and cannot be reused, so the worker
                        // must be replaced before anything else will work. If the GPU was in
                        // use it is the likely culprit (a lost device makes every readback
                        // fail), so retry the same file once on the CPU rather than losing
                        // the run.
                        dispatch("recoverFromEsmPanic", {
                            sampleName: d.sampleName,
                            wasGpu: d.wasGpu,
                        });
                    } else {
                        commit("setEsmError", d.message ?? "generic");
                    }
                } else if (d.gpuEvent) {
                    // Rust has no portable clock, so the gap is the time that batch took.
                    const now = performance.now();
                    const delta = esmLastGpuEventAt > 0
                        ? ` [+${Math.round(now - esmLastGpuEventAt)} ms]`
                        : "";
                    esmLastGpuEventAt = now;
                    console.warn(`[Sparrowhawk] GPU: ${d.message}${delta}`);
                    commit("setEsmGpuEvent", d.message ?? "");
                } else if (d.modelLoading) {
                    // A fresh backend has proven nothing, so guard its first run again.
                    esmGpuProven = false;
                    commit("setLoadingEsmModel", d.stage ?? "");
                } else if (d.backendFallback) {
                    commit("setEsmBackendFallback", d.reason ?? "unknown");
                } else if (d.modelLoaded) {
                    commit("setEsmModelLoaded", {
                        fileName: d.fileName,
                        info: d.info,
                        backend: d.backend,
                    });
                } else if (d.embedProgress) {
                    esmLastProgressAt = performance.now();
                    commit("setEmbeddingProgress", {done: d.done, total: d.total});
                } else if (d.embedded) {
                    clearEsmWatchdog();
                    if (d.backend === "webgpu") esmGpuProven = true;
                    commit("removeEmbeddingFile", d.sampleName);
                    commit("saveEmbeddingResult", {
                        sampleName: d.sampleName,
                        ids: d.ids,
                        lengths: d.lengths,
                        truncated: d.truncated,
                        dim: d.dim,
                        nSequences: d.nSequences,
                        backend: d.backend,
                        elapsedMs: d.elapsedMs,
                        batchMin: d.batchMin,
                        batchMax: d.batchMax,
                        vectors: d.vectors,
                        coords: d.coords,
                    });
                }
            };
        });
        console.log(`Spawned ${workerCount} embedder worker(s)`);
        commit("SET_WORKERS_ESM", pool);
    },

    
    async embedFile(context: ActionContext<RootState, RootState>, payload: {
        acceptFiles: Array<File>,
        use_gpu: boolean,
        gpu_power_pref: number,
        gpu_tasks_max: number,
    }) {
        const {commit, dispatch, state} = context;
        commit("resetAllResults_esm");

        const pool = state.workerState.workers_esm;
        if (!pool.length) {
            console.log("No embedder worker available");
            commit("setEsmError", "worker");
            return;
        }
        if (payload.acceptFiles.length !== 1) {
            commit("setEsmError", "file_count");
            return;
        }

        const file = payload.acceptFiles[0];
        if (!regExpForAnyProteinFasta.test(file.name)) {
            commit("setEsmError", "fasta");
            return;
        }

        const sampleName = file.name.replace(regExpForAnyProteinFasta, "");
        console.log(`Queuing embedding for sample: ${sampleName}`);
        // Kept so the run can be re-issued on the CPU if the wasm module aborts mid-flight.
        commit("setEsmRetry", {file, sampleName});
        commit("addEmbeddingFile", sampleName);
        pool[0].postMessage({
            embed: true,
            file,
            sampleName,
            use_gpu: payload.use_gpu,
            gpu_power_pref: payload.gpu_power_pref,
            gpu_tasks_max: payload.gpu_tasks_max,
        });

        // Only the first GPU run: a broken pipeline never rejects, it just goes quiet.
        if (payload.use_gpu && !esmGpuProven) {
            esmLastProgressAt = performance.now();
            esmStallTimer = setInterval(() => {
                if (performance.now() - esmLastProgressAt < ESM_STALL_MS) return;
                clearEsmWatchdog();
                console.warn(`[Sparrowhawk] GPU made no progress for ${ESM_STALL_MS / 1000}s. Replacing the worker.`);
                commit("removeEmbeddingFile", sampleName);
                dispatch("recoverFromEsmPanic", {sampleName, wasGpu: true});
            }, 2000);
        }
    },


    async recoverFromEsmPanic(context: ActionContext<RootState, RootState>, payload: {
        sampleName?: string,
        wasGpu?: boolean,
    }) {
        const {commit, dispatch, state} = context;
        clearEsmWatchdog();

        // The event channel carries backend-agnostic progress too, so only blame the GPU
        // when the run was actually on it.
        const cause = state.allResults_esm.gpuEvent
            ? `${payload.wasGpu ? "the GPU" : "the engine"} reported: ${state.allResults_esm.gpuEvent}`
            : payload.wasGpu
                ? "the GPU driver dropped the device"
                : "the engine stopped without a reason";
        console.warn(`[Sparrowhawk] embedder wasm aborted; ${cause}. Replacing the worker.`);

        // Take the retry context before anything else, and clear it: whatever happens next
        // must not be able to trigger a second retry.
        const retry = {...state.esmRetry};
        commit("clearEsmRetry");

        // Fresh worker: the old module is poisoned and every call into it would fail.
        await dispatch("initEmbedderWorkers", 1);

        if (!payload.wasGpu || !retry.file) {
            commit("setEsmError", payload.wasGpu ? "gpu_lost" : "wasm_panic");
            return;
        }

        commit("setEsmBackendFallback", cause);
        commit("addEmbeddingFile", retry.sampleName);
        state.workerState.workers_esm[0].postMessage({
            embed: true,
            file: retry.file,
            sampleName: retry.sampleName,
            // The point of the retry: no GPU this time.
            use_gpu: false,
            gpu_power_pref: 0,
            gpu_tasks_max: 0,
        });
    },

    async listGpuAdapters(context: ActionContext<RootState, RootState>) {
        context.commit("setGpuAdapters", await probeGpuAdapters());
    },

    async resetAllResults_esm(context: ActionContext<RootState, RootState>) {
        context.commit("resetAllResults_esm");
    },
};
