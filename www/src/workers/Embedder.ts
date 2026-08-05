import { gunzipSync } from "fflate";
import type { EsmEmbedder } from "@/pkg_esm";
import { ESM_WEIGHTS_FILE_NAME, loadEsmEncoder, loadEsmWeights } from "@/workers/esmWeights";

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

export interface EmbedderOptions {
    use_gpu: boolean;
    gpu_power_pref: number;
    /** 0 = the built-in default; see GPU_TASKS_MAX in the Rust crate. */
    gpu_tasks_max: number;
}

// A Rust panic surfaces as `RuntimeError: unreachable executed` and poisons the module, so
// the worker has to be replaced rather than retried.
// TODO: add to other tabs as well, to catch errors!!
function isWasmPanic(error: unknown): boolean {
    if (error instanceof WebAssembly.RuntimeError) return true;
    const message = error instanceof Error ? error.message : String(error);
    return /unreachable|RuntimeError|recursive use of an object/i.test(message);
}


export class Embedder {
    worker: Worker;
    wasm: WasmModuleAny | null = null;
    embedder: EsmEmbedder | null = null;
    wasmPromise: Promise<WasmModuleAny>;
    private embedderPromise: Promise<void> | null = null;
    private activeBackend: string | null = null;
    private loadedWithGpu: boolean | null = null;
    private currentSample: string | null = null;

    isUsingGpu(): boolean {
        return this.activeBackend === "webgpu";
    }

    getCurrentSample(): string | null {
        return this.currentSample;
    }

    constructor(worker: Worker) {
        this.worker = worker;
        this.wasmPromise = new Promise((resolve, reject) => {
            import("@/pkg_esm")
                .then((w) => {
                    this.wasm = w;
                    if (this.wasm.init_panic_hook) {
                        this.wasm.init_panic_hook();
                    }
                    resolve(w);
                })
                // Without this the rejection reaches `unhandledrejection`, which reports a
                // failed download as a wasm panic and respawns a worker that fails the same way.
                .catch((cause) => {
                    reject(new Error(
                        "Could not load the embedding engine (11 MB wasm module): " +
                        (cause instanceof Error ? cause.message : String(cause)),
                    ));
                });
        });

        // This avoids a failed download (of the wasm binary) to be seen as a wasm (execution) error, by preserving a premise
        // itself that later can be catched in the embedding function (below)
        this.wasmPromise.catch(() => { /* reported when awaited */ });
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    // Like for the AMR detector
    private async ensureEmbedder(useGpu: boolean, powerPref: number, tasksMax: number): Promise<void> {
        if (this.embedder !== null && this.loadedWithGpu === useGpu) return;
        if (this.embedderPromise !== null) return this.embedderPromise;

        this.embedderPromise = this.loadEmbedder(useGpu, powerPref, tasksMax);
        try {
            await this.embedderPromise;
        } finally {
            this.embedderPromise = null;
        }
    }


    private async loadEmbedder(useGpu: boolean, powerPref: number, tasksMax: number): Promise<void> {
        this.worker.postMessage({ modelLoading: true, stage: "engine" });
        const wasm = await this.waitForWasm();

        this.worker.postMessage({ modelLoading: true, stage: "download" });
        const weights = await loadEsmWeights();
        const encoder = await loadEsmEncoder();
        this.worker.postMessage({ modelLoading: true, stage: "init" });

        if (useGpu) {
            try {
                // Also compiles shaders, so this takes a second or two.
                this.embedder = await wasm.EsmEmbedder.create(weights, encoder, true, powerPref, tasksMax);
                this.activeBackend = "webgpu";
            } catch (error) {
                console.warn("[Sparrowhawk] WebGPU init failed, falling back to CPU:", error);
                this.embedder = await wasm.EsmEmbedder.create(weights, encoder, false, 0, 0);
                this.activeBackend = "cpu";
                this.worker.postMessage({
                    backendFallback: true,
                    reason: error instanceof Error ? error.message : String(error),
                });
            }
        } else {
            this.embedder = await wasm.EsmEmbedder.create(weights, encoder, false, 0, 0);
            this.activeBackend = "cpu";
        }

        this.loadedWithGpu = useGpu;
        this.worker.postMessage({
            modelLoaded: true,
            fileName: ESM_WEIGHTS_FILE_NAME,
            backend: this.activeBackend,
            // Assigned in every branch above; TS cannot narrow a mutable field across them.
            info: this.embedder!.info(),
        });
    }

    async embedThisFile(file: File, sampleName: string, opts: EmbedderOptions): Promise<void> {
        this.currentSample = sampleName;
        try {
            await this.ensureEmbedder(opts.use_gpu, opts.gpu_power_pref, opts.gpu_tasks_max);

            // read fasta
            const raw = new Uint8Array(await file.arrayBuffer());
            const fastaBytes = file.name.endsWith(".gz") ? gunzipSync(raw) : raw;

            const t0 = performance.now();
            const batch = await this.embedder!.embed_fasta(fastaBytes);

            const meta = JSON.parse(batch.meta);
            // The getters copy out of wasm memory, so these are safe to transfer below.
            const vectors: Float32Array = batch.vectors;
            const coords: Float32Array = batch.coords;
            // Without this every file leaks n * 320 floats in the wasm heap.
            batch.free();

            this.worker.postMessage(
                {
                    embedded: true,
                    sampleName,
                    backend: this.activeBackend,
                    ids: meta.ids,
                    lengths: meta.lengths,
                    truncated: meta.truncated,
                    dim: meta.dim,
                    nSequences: meta.n,
                    elapsedMs: Math.round(performance.now() - t0),
                    batchMin: meta.batch_min,
                    batchMax: meta.batch_max,
                    vectors,
                    coords,
                },
                // Transfer, not clone: n * 320 * 4 bytes.
                [vectors.buffer, coords.buffer],
            );
        } catch (error) {
            if (isWasmPanic(error)) {
                this.worker.postMessage({
                    error: true,
                    wasmPanic: true,
                    // Only a GPU run loses the device, so say whether a CPU retry is worth it.
                    wasGpu: this.activeBackend === "webgpu",
                    sampleName,
                    message: error instanceof Error ? error.message : String(error),
                });
                return;
            }
            this.worker.postMessage({
                error: true,
                sampleName,
                message: error instanceof Error ? error.message : String(error),
            });
        }
    }

    // We don't erase the model to re-download it!
    resetAll(): void {
        /* keep the model resident */
    }
}
