import { loadWasmPackage } from "@/wasm/loader";
import { skaPackage } from "@/wasm/packages/ska";
import type { TreeResult } from "@/types";
import type { AlignData, SkaData } from "@/pkg_ska";

// These interfaces are those of the JSONs for retrieving info and returning info
interface MapResult {
    "Number of variants": number;
    "Coverage": number;
    "Mapped sequences": string[];
    "VCF": string;
}

interface AlignResult extends TreeResult {
    alignmentAvailable: boolean;
}

interface ClusterLabels {
    [sampleName: string]: number;
}

interface TransmissionGraphNode { id: string; cluster: number; }
interface TransmissionGraphLink { source: string; target: string; snp_distance: number; }
interface TransmissionGraphData {
    nodes: TransmissionGraphNode[];
    links: TransmissionGraphLink[];
}

type SkaModule = typeof import("@/pkg_ska");

export class Mapper {
    worker: Worker;
    wasm: SkaModule | null;
    SkaData: SkaData | null;
    AlignData: AlignData | null;
    wasmPromise: Promise<SkaModule>;
    wasmMemory: WebAssembly.Memory | null = null;
    private activeAlignmentRunId: number | null = null;
    private alignmentStartAt = 0;
    private completedAlignmentRunId: number | null = null;
    private alignmentK: number | null = null;
    private alignmentRc = false;
    private alignmentFrozen = false;
    private alignmentExportPartSize: number | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
        this.SkaData = null;
        this.AlignData = null;
        this.wasm = null;
        // We use now loadWasmPackage to be able to get either the wasm64 or wasm32 one.
        this.wasmPromise = loadWasmPackage(skaPackage, (notice) => worker.postMessage(notice))
            .then((loaded) => {
                this.wasm = loaded.api;
                this.wasmMemory = loaded.memory ?? null;
                return loaded.api;
            });
        // The caller reports loading failures when it awaits this promise.
        void this.wasmPromise.catch(() => { /* handled by the requesting operation */ });
    }

    waitForWasm(): Promise<SkaModule> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async set_ref(file: File, k: number, rc: boolean, ambig_mask: boolean, repeat_mask: boolean): Promise<void> {
        const wasm = await this.waitForWasm();

        if (this.SkaData === null) {
            this.SkaData = wasm.SkaData.new(file, k, rc, ambig_mask, repeat_mask);
        }
        this.worker.postMessage({
            ref: file,
            sequences: this.SkaData!.get_reference().split('\n')
        });
    }

    map(file: File, revReadFile: File | null, proportion_reads: number, min_count: number, min_qual: number, qual_filter: number): void {
        console.log("Mapping reads to reference with proportion_reads: " + proportion_reads);
        if (this.SkaData === null) {
            throw new Error("SkaRef::map - reference does not exist yet.");
        }

        try {
            const t0 = performance.now();
            const outname = (revReadFile != null) ? file.name.replace(new RegExp("(?:_1)?\\.(?:fa|fna|fasta|fq|fnq|fastq)(?:\\.gz)?" + String.fromCharCode(36)), "") : file.name;
            const results: MapResult = JSON.parse(this.SkaData.map(file, revReadFile, proportion_reads, min_count, min_qual, qual_filter, outname));

            this.worker.postMessage({
                nb_variants: results["Number of variants"],
                coverage: results["Coverage"],
                name: outname,
                mapped_sequences: results["Mapped sequences"],
                mapping_vcf: results["VCF"],
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes(),
            });
        } catch {
            this.worker.postMessage({ error: true, message: 'memory' });
        }
    }

    async beginAlignment(runId: number, fileNames: string[], k: number, rc: boolean, append: boolean): Promise<void> {
        this.activeAlignmentRunId = runId;
        let candidate: AlignData | null = null;
        let planAccepted = false;
        try {
            const wasm = await this.waitForWasm();
            if (this.activeAlignmentRunId !== runId) return;
            if (append) {
                if (!this.AlignData || this.completedAlignmentRunId === null) {
                    throw new Error("No completed dataset is available for additions. Clear results and start again.");
                }
                if (this.alignmentFrozen) throw new Error("The alignment is closed to new samples.");
                if (this.alignmentK !== k || this.alignmentRc !== rc) {
                    throw new Error("k and reverse-complement settings must match the existing dataset.");
                }
                candidate = this.AlignData;
            } else {
                candidate = wasm.AlignData.new(k, runId);
            }
            const samples = JSON.parse(candidate!.plan_samples(fileNames, runId));
            planAccepted = true;
            if (!append) {
                this.AlignData?.free();
                this.AlignData = candidate;
                this.alignmentFrozen = false;
                this.alignmentK = k;
                this.alignmentRc = rc;
            }
            this.alignmentExportPartSize = null;
            this.completedAlignmentRunId = null;
            this.alignmentStartAt = performance.now();
            this.worker.postMessage({ alignmentPlan: true, runId, samples });
        } catch (error) {
            if (candidate && candidate !== this.AlignData) candidate.free();
            if (this.activeAlignmentRunId === runId) {
                if (!planAccepted) this.activeAlignmentRunId = null;
                this.reportAlignmentError(runId, error, !planAccepted);
            }
        }
    }

    addAlignmentSample(
        runId: number,
        sampleIndex: number,
        sampleName: string,
        packedKmers: Uint32Array,
    ): void {
        if (this.activeAlignmentRunId !== runId || this.AlignData === null) return;

        try {
            this.AlignData.add_sample_words(sampleIndex, sampleName, packedKmers);
            this.worker.postMessage({ alignmentSampleAdded: true, runId, sampleIndex });
        } catch (error) {
            this.reportAlignmentError(runId, error);
        }
    }

    finishAlignment(runId: number): void {
        const alignData = this.AlignData;
        if (this.activeAlignmentRunId !== runId || alignData === null) return;

        try {
            const results: AlignResult = JSON.parse(alignData.finish_alignment());
            const distancesCsvGzip = alignData.get_distances_csv_gzip();
            const distancesCsvGzipBytes = distancesCsvGzip.byteLength;
            const elapsedMs = Math.round(performance.now() - this.alignmentStartAt);
            this.worker.postMessage({
                aligned: true,
                runId,
                names: results.names,
                treeTipNames: results.treeTipNames,
                newick: results.newick,
                alignmentAvailable: results.alignmentAvailable,
                alignmentFrozen: false,
                k: this.alignmentK,
                rc: this.alignmentRc,
                distances_csv_gzip: distancesCsvGzip,
                elapsedMs,
                wasmMemoryBytes: this.memoryBytes(),
            }, [distancesCsvGzip.buffer as ArrayBuffer]);
            console.info("[ska-align] Alignment completed", {
                runId,
                sampleCount: results.names.length,
                distancesCsvGzipBytes,
                elapsedMs,
                wasmMemoryBytes: this.memoryBytes(),
            });
            this.completedAlignmentRunId = runId;
            this.activeAlignmentRunId = null;
        } catch (error) {
            this.reportAlignmentError(runId, error);
        }
    }

    beginAlignmentExport(runId: number, requestId: number, samplesPerPart: number): void {
        try {
            if (this.completedAlignmentRunId !== runId || !this.AlignData) {
                throw new Error("The requested alignment dataset is no longer available.");
            }
            if (!Number.isInteger(samplesPerPart) || samplesPerPart < 1) {
                throw new Error("The requested alignment export part size is invalid.");
            }
            this.AlignData.begin_alignment_export(samplesPerPart);
            this.alignmentFrozen = true;
            this.alignmentExportPartSize = samplesPerPart;
            console.info("[ska-align] Started batched alignment export", { runId, requestId, samplesPerPart });
            this.nextAlignmentExportChunk(runId, requestId);
        } catch (error) {
            this.reportAlignmentExportError(runId, requestId, error);
        }
    }

    nextAlignmentExportChunk(runId: number, requestId: number): void {
        try {
            if (this.completedAlignmentRunId !== runId || !this.AlignData || this.alignmentExportPartSize === null) {
                throw new Error("No resumable alignment export is available.");
            }
            const alignData = this.AlignData;
            const bytes = alignData.next_alignment_export_chunk();
            const samplesProcessed = alignData.alignment_export_samples_processed();
            const partComplete = alignData.alignment_export_part_complete();
            const complete = alignData.alignment_export_complete();
            this.worker.postMessage({
                alignmentDownloadChunk: true,
                runId,
                requestId,
                bytes,
                samplesProcessed,
                partComplete,
                complete,
                wasmMemoryBytes: this.memoryBytes(),
            }, [bytes.buffer as ArrayBuffer]);
            if (complete) {
                this.alignmentExportPartSize = null;
                console.info("[ska-align] Alignment output generated", { runId, requestId, samplesProcessed, wasmMemoryBytes: this.memoryBytes() });
            }
        } catch (error) {
            this.reportAlignmentExportError(runId, requestId, error);
        }
    }

    cancelAlignmentExport(runId: number, requestId: number): void {
        if (this.completedAlignmentRunId !== runId || !this.AlignData) return;
        this.AlignData.cancel_alignment_export();
        this.alignmentExportPartSize = null;
        this.worker.postMessage({ alignmentExportCancelled: true, runId, requestId });
    }

    finishAlignmentExport(runId: number, requestId: number): void {
        try {
            if (this.completedAlignmentRunId !== runId || !this.AlignData) {
                throw new Error("The requested alignment dataset is no longer available.");
            }
            this.AlignData.finish_alignment_export();
            this.alignmentExportPartSize = null;
            this.worker.postMessage({ alignmentExportFinished: true, runId, requestId });
        } catch (error) {
            this.reportAlignmentExportError(runId, requestId, error);
        }
    }

    private reportAlignmentExportError(runId: number, requestId: number, error: unknown): void {
        const message = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment export failed", { runId, requestId, error });
        this.worker.postMessage({ alignmentDownloadError: true, runId, requestId, message });
    }

    private resetAlignmentExport(): void {
        if (this.AlignData && this.alignmentExportPartSize !== null) {
            this.AlignData.cancel_alignment_export();
        }
        this.alignmentExportPartSize = null;
    }

    cancelAlignment(runId: number): void {
        if (this.activeAlignmentRunId === runId) {
            this.AlignData?.free();
            this.AlignData = null;
            this.completedAlignmentRunId = null;
            this.activeAlignmentRunId = null;
            this.alignmentExportPartSize = null;
        }
    }

    private reportAlignmentError(runId: number, error: unknown, validationOnly = false): void {
        const message = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment failed", { runId, wasmMemoryBytes: this.memoryBytes(), error });
        this.worker.postMessage({ error: true, alignmentError: true, validationOnly, runId, message });
    }

    async clusterUploadedAlignment(file: File, snp_threshold: number, requestId: number): Promise<void> {
        console.log("Running standalone transmission clustering with SNP threshold: " + snp_threshold);
        let alignData: AlignData | null = null;
        try {
            const t0 = performance.now();
            const wasm = await this.waitForWasm();
            this.worker.postMessage({ standalone: true, requestId, stage: "distances" });
            const alignmentText = await file.text();
            const imported: AlignData = wasm.AlignData.from_alignment_text(alignmentText);
            alignData = imported;
            this.worker.postMessage({ standalone: true, requestId, stage: "clustering" });
            const clusters: ClusterLabels = JSON.parse(wasm.ska_cluster(imported, snp_threshold));
            const graph: TransmissionGraphData = JSON.parse(imported.get_graph_json(snp_threshold));
            this.worker.postMessage({ standalone: true, requestId, stage: "tree" });
            let tree: TreeResult | null = null;
            let treeError: string | null = null;
            try {
                tree = JSON.parse(imported.build_tree());
            } catch (error) {
                treeError = error instanceof Error ? error.message : String(error);
            }
            this.worker.postMessage({ clustered: true, standalone: true, requestId, clusters, graph, tree, treeError,
                distancesCapped: imported.distances_capped(),
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes() });
        } catch (error) {
            const message = error instanceof Error ? error.message : String(error);
            this.worker.postMessage({ error: true, clusterError: true, standalone: true, requestId, message });
        } finally {
            alignData?.free();
        }
    }

    cluster(snp_threshold: number, runId: number, requestId: number): void {
        console.log("Running transmission clustering with SNP threshold: " + snp_threshold);
        try {
            if (this.wasm === null || this.AlignData === null || this.completedAlignmentRunId !== runId) throw new Error("No completed alignment dataset is available.");
            const clusters: ClusterLabels = JSON.parse(this.wasm.ska_cluster(this.AlignData, snp_threshold));
            const graph: TransmissionGraphData = JSON.parse(this.AlignData.get_graph_json(snp_threshold));
            this.worker.postMessage({ clustered: true, runId, requestId, clusters, graph });
        } catch (error) {
            this.worker.postMessage({ error: true, clusterError: true, runId, requestId, message: error instanceof Error ? error.message : String(error) });
        }
    }

    resetAll(): void {
        this.SkaData = null;
        this.AlignData?.free();
        this.AlignData = null;
        this.activeAlignmentRunId = null;
        this.completedAlignmentRunId = null;
        this.alignmentFrozen = false;
        this.resetAlignmentExport();
    }
}
