import {ActionContext} from "vuex";
import {RootState} from "@/store/state";
import WorkerSketcher from '@/workers/Sketcher.worker';
import WorkerCaller from '@/workers/Caller.worker';
import {regExpWithTwoNumbers, regExpForAnyFastx, getFilesToProcess} from "@/utils";

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
    }) {
        const {commit, state} = context;
        console.log("Processing query of uploaded files for alignment...")

        // Initialize the aligned state so that we can know that it is loading
        commit("setAligned", {aligned: false, names: [], newick: ""})

        // Set aligning state
        commit("setAligningState", true);

        const messageData = {
            align: true,
            files: payload.acceptFiles,
            k: payload.k,
            proportion_reads: payload.proportion_reads,
            rc: payload.rc,
            min_count: payload.min_count,
            min_qual: payload.min_qual,
            qual_filter: payload.qual_filter,
        };

        if (state.workerState.worker_ska) {
            state.workerState.worker_ska.postMessage(messageData);
            state.workerState.worker_ska.onmessage = (message) => {
                // Clear aligning state
                commit("setAligningState", false);

                if (message.data.error) {
                    commit("setSkaAlignError", message.data.message ?? "generic");
                    return;
                }
                commit("setAligned", message.data);
            };
        }
    },

    async processCluster(context: ActionContext<RootState, RootState>, payload: { snp_threshold: number }) {
        const {commit, state} = context;
        console.log("Running transmission clustering with SNP threshold: " + payload.snp_threshold);

        commit("setClusteringState", true);

        if (state.workerState.worker_ska) {
            state.workerState.worker_ska.postMessage({ cluster: true, snp_threshold: payload.snp_threshold });
            state.workerState.worker_ska.onmessage = (message) => {
                commit("setClusteringState", false);
                if (message.data.error) {
                    commit("setSkaClusterError", message.data.message ?? "generic");
                    return;
                }
                commit("setClusterResults", { clusters: message.data.clusters, graph: message.data.graph });
            };
        }
    },

    async resetAllResults_ska(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
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

        // Attach a result handler to each worker in the pool
        pool.forEach((worker) => {
            worker.onmessage = (messageData) => {
                if (messageData.data instanceof Object) {
                    if (messageData.data.error) {
                        const sampleName = messageData.data.sampleName;
                        commit("setSketchlibError", messageData.data.message ?? "generic");
                        if (sampleName) commit("removeIdentifyingFile", sampleName);
                    } else if ("probs" in messageData.data && "sampleName" in messageData.data) {
                        const sampleName = messageData.data.sampleName;
                        console.log("Saving results for sample: " + sampleName);
                        commit("removeIdentifyingFile", sampleName);
                        commit("saveIDResults", {
                            sampleName,
                            probs: messageData.data.probs,
                            names: messageData.data.names,
                            metadata: messageData.data.metadata,
                        });
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
        // Terminate existing workers
        for (const worker of state.workerState.workers_sketchlib) {
            worker.terminate();
        }
        // Spawn new pool
        const pool: Worker[] = [];
        for (let i = 0; i < numWorkers; i++) {
            pool.push(new WorkerSketcher());
        }
        console.log(`Spawned ${numWorkers} sketchlib worker(s)`);
        commit("SET_WORKERS_SKETCHLIB", pool);
    },

    async resetAllResults_sketchlib(context: ActionContext<RootState, RootState>) {
        const {commit} = context;
        commit("resetAllResults_sketchlib");
    },

    // ORPHOS
    async initCallerWorkers(context: ActionContext<RootState, RootState>, numWorkers: number) {
        const { commit, state } = context;
        for (const worker of state.workerState.workers_orphos) {
            worker.terminate();
        }
        const pool: Worker[] = [];
        for (let i = 0; i < numWorkers; i++) {
            pool.push(new WorkerCaller());
        }
        pool.forEach(worker => {
            worker.onmessage = (msg) => {
                if (msg.data?.error) {
                    commit("setOrphosError", msg.data.message ?? "generic");
                    if (msg.data.fileName) commit("removeCallingGenesFile", msg.data.fileName);
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
                    });
                }
            };
        });
        console.log(`Spawned ${numWorkers} caller worker(s)`);
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
        }) {
        const { commit, state } = context;
        const pool = state.workerState.workers_orphos;
        if (!pool.length) {
            console.log("No caller workers available");
            return;
        }
        const indxlist = getFilesToProcess(payload.acceptFiles);
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
    async loadAmrIndex(context: ActionContext<RootState, RootState>, payload: { file: File }) {
        const { commit, state } = context;
        if (!state.workerState.worker_amr) return;
        commit("setLoadingAmrIndex");
        state.workerState.worker_amr.postMessage({ loadIndex: true, file: payload.file });
        state.workerState.worker_amr.onmessage = (msg) => {
            if (msg.data.error) {
                state.allResults_amr.isLoadingIndex = false;
                commit("setAmrError", msg.data.message ?? "index");
            } else if (msg.data.indexLoaded) {
                commit("setAmrIndexLoaded", { fileName: msg.data.fileName, info: msg.data.info });
            }
        };
    },

    async detectAmr(context: ActionContext<RootState, RootState>, payload: { files: Array<File>; min_gene_hits: number; min_family_hits: number }) {
        const { commit, state } = context;
        if (!state.workerState.worker_amr) return;
        const workerAmr = state.workerState.worker_amr;
        workerAmr.onmessage = (msg) => {
            if (msg.data.error) {
                commit("setAmrError", msg.data.message ?? "generic");
                if (msg.data.sampleName) commit("removeDetectingAmrFile", msg.data.sampleName);
            } else if (msg.data.detected) {
                commit("removeDetectingAmrFile", msg.data.sampleName);
                commit("saveAmrResult", msg.data.result);
            }
        };

        for (const file of payload.files) {
            const sampleName = file.name.replace(regExpForAnyFastx, "");
            commit("addDetectingAmrFile", sampleName);
            workerAmr.postMessage({
                detect: true,
                file,
                sampleName,
                min_gene_hits: payload.min_gene_hits,
                min_family_hits: payload.min_family_hits,
            });
        }
    },

    async resetAllResults_amr(context: ActionContext<RootState, RootState>) {
        context.commit("resetAllResults_amr");
    },
};
