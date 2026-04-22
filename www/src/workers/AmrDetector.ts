import type { WasmAmrIndex } from "@/pkg_amr";

export class AmrDetector {
    worker: Worker;
    index: WasmAmrIndex | null = null;
    wasm: typeof import("@/pkg_amr") | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
    }

    async loadIndex(file: File): Promise<void> {
        try {
            const wasm = await import("@/pkg_amr");
            this.wasm = wasm;
            const buf = await file.arrayBuffer();
            this.index = new wasm.WasmAmrIndex(new Uint8Array(buf));
            this.worker.postMessage({ indexLoaded: true, fileName: file.name, info: this.index.info() });
        } catch {
            this.worker.postMessage({ error: true, message: "index" });
        }
    }

    async detect(file: File, sampleName: string, min_gene_hits: number, min_family_hits: number): Promise<void> {
        try {
            if (!this.index || !this.wasm) {
                throw new Error("AMR index is not loaded");
            }
            const buf = new Uint8Array(await file.arrayBuffer());
            const session = new this.wasm.WasmAmrSession(this.index, min_gene_hits, min_family_hits);
            session.push_chunk(buf);
            const result = session.finish(sampleName) as any;
            this.worker.postMessage({ detected: true, sampleName, result });
        } catch {
            this.worker.postMessage({ error: true, sampleName, message: "memory" });
        }
    }

    resetAll(): void {
        this.index = null;
        this.wasm = null;
    }
}
