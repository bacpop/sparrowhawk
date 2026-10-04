import {loadAssetBlob} from "@/platform/files";
import type { SketchlibData } from "@/pkg_sketchlib";

interface IdentifyResult {
    ani: number[];
    ranks: number[];
    names: string[];
    metadata: string[];
}

type SketchlibModule = typeof import("@/pkg_sketchlib");

export class Sketcher {
    worker: Worker;
    wasm: SketchlibModule | null;
    SketchlibData: SketchlibData | null;
    wasmPromise: Promise<SketchlibModule>;
    wasmMemory: WebAssembly.Memory | null = null;

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
        import("@/pkg_sketchlib/index_bg.wasm").then((m) => { this.wasmMemory = m.memory; });
    }

    waitForWasm(): Promise<SketchlibModule> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async identifyThisFile(file1: File, file2: File | null, sampleName: string, proportion_reads: number, min_count: number, min_qual: number): Promise<void> {
        console.log("Starting identification for sample: " + sampleName);
        const wasm = await this.waitForWasm();

        try {
            if (this.SketchlibData === null) {
                let invertedindex: Blob;
                try {
                    invertedindex = await loadAssetBlob("inverted_k_17_ss_50.ski");
                } catch (error) {
                    console.error("Failed to load sketchlib asset", error);
                    this.worker.postMessage({error: true, sampleName, message: "asset"});
                    return;
                }

                // Typed as File, but the Rust reader only uses size and slice, which a Blob has too.
                this.SketchlibData = wasm.SketchlibData.new(invertedindex as File);
            }

            this.SketchlibData!.query(file1, file2, proportion_reads, min_count, min_qual);

            const results: IdentifyResult = JSON.parse(this.SketchlibData!.get_ani(3));
            const ani = results.ani;

            this.worker.postMessage({
                sampleName: sampleName,
                ani,
                ranks: results.ranks,
                names: results.names,
                metadata: results.metadata,
                wasmMemoryBytes: this.memoryBytes()
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
