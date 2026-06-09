interface IdentifyResult {
    ani: number[];
    ranks: number[];
    names: string[];
    metadata: string[];
}

interface SketchlibData {
    query(file1: File, file2: File | null, proportion_reads: number, min_count: number, min_qual: number): Promise<void>;
    get_ani(top_n: number): string;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

export class Sketcher {
    worker: Worker;
    wasm: WasmModuleAny | null;
    SketchlibData: SketchlibData | null;
    wasmPromise: Promise<WasmModuleAny>;

    constructor(worker: Worker) {
        this.worker = worker;
        this.SketchlibData = null;
        this.wasm = null;
        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg_sketchlib")
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

    async identifyThisFile(file1: File, file2: File | null, sampleName: string, proportion_reads: number, min_count: number, min_qual: number): Promise<void> {
        console.log("Starting identification for sample: " + sampleName);
        await this.waitForWasm();

        try {
            if (this.SketchlibData === null) {
                const response = await fetch('/inverted_k_17_ss_50.ski');
                const invertedindex = await response.blob();

                this.SketchlibData = await this.wasm.SketchlibData.new(invertedindex);
            }

            await this.SketchlibData!.query(file1, file2, proportion_reads, min_count, min_qual);

            const results: IdentifyResult = JSON.parse(this.SketchlibData!.get_ani(3));
            const ani = results.ani;

            this.worker.postMessage({
                sampleName: sampleName,
                ani,
                ranks: results.ranks,
                names: results.names,
                metadata: results.metadata
            });
        } catch (error) {
            this.worker.postMessage({
                error: true,
                sampleName,
                message: error instanceof Error ? error.message : String(error),
            });
        }
    }

    resetAll(): void {
        this.SketchlibData = null;
    }
}
