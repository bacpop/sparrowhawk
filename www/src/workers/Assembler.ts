import type { AssemblyHelper } from "@/pkg";
import { loadWasmPackage } from "@/wasm/loader";
import { sparrowhawkPackage } from "@/wasm/packages/sparrowhawk";

interface PreprocessingResult {
    nkmers: number;
    histo: number[];
    used_min_count: number;
    k: number;
    multik: boolean;
}

interface AssemblyResult {
    ncontigs: number;
}

type AsmModule = typeof import("@/pkg");

/** What the page is told when a call fails. A WebAssembly trap (a Rust panic, or memory running out)
 *  is a crash; anything else is an error the assembler raised on purpose, whose message says why. */
function failure(error: unknown): { reset: true, error: string, crashed: boolean } {
    return {
        reset: true,
        error: error instanceof Error ? error.message : String(error),
        crashed: error instanceof WebAssembly.RuntimeError,
    };
}

export class Assembler {
    worker: Worker;
    wasm: AsmModule | null;
    helper: AssemblyHelper | null;
    wasmPromise: Promise<AsmModule>;
    wasmMemory: WebAssembly.Memory | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
        this.helper = null;
        this.wasm = null;

        this.wasmPromise = loadWasmPackage(sparrowhawkPackage, (notice) => worker.postMessage(notice))
            .then((loaded) => {
                this.wasm = loaded.api;
                this.wasmMemory = loaded.memory ?? null;
                return loaded.api;
            });
        // The caller reports loading failures when it awaits this promise.
        void this.wasmPromise.catch(() => { /* handled by the requesting operation */ });
    }

    waitForWasm(): Promise<AsmModule> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async preprocess(
        file1: File,
        file2: File | null,
        k: number,
        verbose: boolean,
        min_count: number,
        min_qual: number | null,
        do_bloom: boolean,
        do_fit: boolean,
        multik: boolean,
        ks: number[],
        min_contig_length: number
    ): Promise<void> {
        const t0 = performance.now();

        if (this.helper === null) {
            try {
                const wasm = await this.waitForWasm();
                this.helper = wasm.AssemblyHelper.new(
                    k, verbose,
                    min_count, min_qual,
                    do_bloom, do_fit,
                    min_contig_length
                );
                if (multik) {
                    this.helper.set_multik(new Uint32Array(ks));
                }
                this.helper.preprocess(file1, file2);
            } catch (error) {
                console.error(error);
                this.worker.postMessage(failure(error));
                return;
            }
        }

        const resultsjson: PreprocessingResult = JSON.parse(
            this.helper!.get_preprocessing_info()
        );

        this.worker.postMessage({
            nKmers:         resultsjson.nkmers,
            histo:          resultsjson.histo,
            used_min_count: resultsjson.used_min_count,
            k:              resultsjson.k,
            multik:         resultsjson.multik,
            elapsedMs:      Math.round(performance.now() - t0),
        });
    }

    async assemble(): Promise<void> {
        console.log("Initiating assembly from worker");
        const t0 = performance.now();
        try {
            this.helper!.assemble();
            const resultsjson: AssemblyResult = JSON.parse(this.helper!.get_assembly());
            const fastaGz = this.helper!.take_fasta_gz();
            const dotGz = this.helper!.take_dot_gz();
            const gfaGz = this.helper!.take_gfa_gz();
            const gfa2Gz = this.helper!.take_gfa2_gz();

            // Transfer the buffers: the main thread keeps the compressed outputs for repeat downloads.
            this.worker.postMessage({
                ncontigs: resultsjson.ncontigs,
                fastaGz, dotGz, gfaGz, gfa2Gz,
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes(),
            }, [fastaGz.buffer as ArrayBuffer, dotGz.buffer as ArrayBuffer,
                gfaGz.buffer as ArrayBuffer, gfa2Gz.buffer as ArrayBuffer]);
        } catch (error) {
            console.error(error);
            this.worker.postMessage(failure(error));
        }
    }

    /** Free the helper now: it holds a whole counted k, which would otherwise wait for the GC. */
    resetAll(): void {
        this.helper?.free();
        this.helper = null;
    }
}
