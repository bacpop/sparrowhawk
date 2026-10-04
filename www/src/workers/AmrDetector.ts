import { gunzipSync } from "fflate";
import { AMR_INDEX_FILE_NAME, createAmrDetector } from "@/workers/amrIndex";

import type { AmrDetector } from "@/pkg_amr";

type AmrModule = typeof import("@/pkg_amr");

export class AmrDetectorWorker {
    worker: Worker;
    wasm: AmrModule | null;
    detector: AmrDetector | null;
    wasmPromise: Promise<AmrModule>;
    detectorPromise: Promise<void> | null;

    wasmMemory: WebAssembly.Memory | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
        this.wasm = null;
        this.detector = null;
        this.detectorPromise = null;
        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg_amr")
                .then((w) => {
                    this.wasm = w;
                    if (this.wasm.init_panic_hook) {
                        this.wasm.init_panic_hook();
                    }
                    resolve(w);
                });
        });
        import("@/pkg_amr/index_bg.wasm").then((m) => { this.wasmMemory = m.memory; });
    }

    waitForWasm(): Promise<AmrModule> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async ensureDetector(): Promise<void> {
        if (this.detector !== null) return;
        if (this.detectorPromise !== null) return this.detectorPromise;

        this.detectorPromise = this.loadDetector();
        try {
            await this.detectorPromise;
        } finally {
            this.detectorPromise = null;
        }
    }

    async detectThisFile(file: File, sampleName: string, min_gene_fraction: number, min_gene_group_fraction: number): Promise<void> {
        try {
            await this.ensureDetector();
            const raw = new Uint8Array(await file.arrayBuffer());
            const fastaBytes = file.name.endsWith(".gz") ? gunzipSync(raw) : raw;
            const t0 = performance.now();
            const json = this.detector!.detect_direct(
                sampleName,
                fastaBytes,
                min_gene_fraction,
                min_gene_group_fraction
            );
            this.worker.postMessage({
                detected: true,
                sampleName,
                result: { ...JSON.parse(json), elapsedMs: Math.round(performance.now() - t0), wasmMemoryBytes: this.memoryBytes() },
            });
        } catch (error) {
            this.worker.postMessage({
                error: true,
                sampleName,
                message: error instanceof Error ? error.message : String(error),
            });
        }
    }

    private async loadDetector(): Promise<void> {
        const wasm = await this.waitForWasm();
        this.detector = await createAmrDetector(wasm);
        this.worker.postMessage({
            indexLoaded: true,
            fileName: AMR_INDEX_FILE_NAME,
            info: this.detector.info(),
        });
    }
}
