import type { WasmIndex } from "@/pkg_deacon";

interface FilterStats {
    readsIn: number;
    readsOut: number;
    basesIn: number;
    basesOut: number;
}

function concat(chunks: Uint8Array[]): Uint8Array {
    const total = chunks.reduce((n, c) => n + c.length, 0);
    const out = new Uint8Array(total);
    let offset = 0;
    for (const c of chunks) {
        out.set(c, offset);
        offset += c.length;
    }
    return out;
}

export class Depleter {
    worker: Worker;
    index: WasmIndex | null = null;
    wasm: typeof import("@/pkg_deacon") | null = null;
    wasmMemory: WebAssembly.Memory | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async loadIndex(file: File): Promise<void> {
        try {
            const wasm = await import("@/pkg_deacon");
            this.wasm = wasm;
            import("@/pkg_deacon/index_bg.wasm").then((m) => { this.wasmMemory = m.memory; });
            const buf = await file.arrayBuffer();
            this.index = new wasm.WasmIndex(new Uint8Array(buf));
            this.worker.postMessage({ indexLoaded: true, fileName: file.name, info: this.index.info() });
        } catch {
            this.worker.postMessage({ error: true, message: 'index' });
        }
    }

    async filterReads(file: File, revFile: File | null, deplete: boolean, abs_threshold: number, rel_threshold: number, sampleName: string): Promise<void> {
        const isGz = (f: File) => f.name.endsWith('.gz');

        const runSession = async (f: File): Promise<{ out: Uint8Array; stats: FilterStats }> => {
            const buf = new Uint8Array(await f.arrayBuffer());
            const session = new this.wasm!.WasmFilterSession(
                this.index!, deplete, abs_threshold, rel_threshold, isGz(f), true
            );
            const chunks: Uint8Array[] = [];
            const part = session.push_chunk(buf);
            if (part.length > 0) chunks.push(part);
            const tail = session.finish();
            if (tail.length > 0) chunks.push(tail);
            const stats: FilterStats = session.stats();
            const out = concat(chunks);
            return { out, stats };
        };

        try {
            const t0 = performance.now();
            const { out: outputGzip, stats: s1 } = await runSession(file);
            let outputGzip2: Uint8Array | null = null;
            let s2stats = { readsIn: 0, readsOut: 0 };

            if (revFile) {
                const { out, stats } = await runSession(revFile);
                outputGzip2 = out;
                s2stats = stats;
            }

            const total   = s1.readsIn + s2stats.readsIn;
            const kept    = s1.readsOut + s2stats.readsOut;
            const removed = total - kept;

            this.worker.postMessage(
                { filtered: true, sampleName, total, kept, removed, outputGzip, outputGzip2,
                  elapsedMs: Math.round(performance.now() - t0),
                  wasmMemoryBytes: this.memoryBytes() },
                outputGzip2
                    ? [outputGzip.buffer, outputGzip2.buffer]
                    : [outputGzip.buffer]
            );
        } catch {
            this.worker.postMessage({ error: true, sampleName, message: 'memory' });
        }
    }

    resetAll(): void {
        this.index = null;
        this.wasm = null;
    }
}
