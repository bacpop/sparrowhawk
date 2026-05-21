import { gunzipSync } from "fflate";

const AMR_INDEX_FILE_NAME = "amrfinderplus_2026-03-24.1_dna_k31.amridx.gz";
const AMR_INDEX_URL = `/${AMR_INDEX_FILE_NAME}`;

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

export class AmrDetectorWorker {
    worker: Worker;
    wasm: WasmModuleAny | null;
    detector: WasmModuleAny | null;
    wasmPromise: Promise<WasmModuleAny>;
    detectorPromise: Promise<void> | null;

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
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
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

    async detectThisFile(file: File, sampleName: string, min_gene_fraction: number, min_family_fraction: number): Promise<void> {
        try {
            await this.ensureDetector();
            const raw = new Uint8Array(await file.arrayBuffer());
            const fastaBytes = file.name.endsWith(".gz") ? gunzipSync(raw) : raw;
            const json = this.detector!.detect_direct(
                sampleName,
                fastaBytes,
                min_gene_fraction,
                min_family_fraction
            );
            this.worker.postMessage({
                detected: true,
                sampleName,
                result: JSON.parse(json),
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
        const response = await fetch(AMR_INDEX_URL);
        if (!response.ok) {
            throw new Error("index");
        }
        const compressed = new Uint8Array(await response.arrayBuffer());
        const indexBytes = gunzipSync(compressed);
        this.detector = new wasm.AmrDetector(indexBytes);
        this.worker.postMessage({
            indexLoaded: true,
            fileName: AMR_INDEX_FILE_NAME,
            info: this.detector.info(),
        });
    }
}
