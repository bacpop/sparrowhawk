import { createAmrDetector } from "@/workers/amrIndex";
import { buildGeneCallingAmrTsv } from "@/amrTsv";
import { AmrDetectionResult, GeneMetadataMap } from "@/types";
interface CallResult {
    output_file: string;
    gene_count: number;
    sequence_count: number;
}

interface OrphosData {
    read_fasta(input_file: File): void;
    index_fasta(): void;
    call_genes(): void;
    get_results(format: string): string;
    get_annotated_results(format: string, amr_json: string): string;
    get_cds_fasta(): string;
    get_gene_metadata_json(): string;
    take_fasta_bgz(): Uint8Array<ArrayBuffer>;
    take_fasta_fai(): Uint8Array<ArrayBuffer>;
    take_fasta_gzi(): Uint8Array<ArrayBuffer>;
    take_gff_bgz():   Uint8Array<ArrayBuffer>;
    take_gff_csi():   Uint8Array<ArrayBuffer>;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

export class Caller {
    worker: Worker;
    wasm: WasmModuleAny | null;
    OrphosData: OrphosData | null;
    wasmPromise: Promise<WasmModuleAny>;
    amrWasm: WasmModuleAny | null;
    amrWasmPromise: Promise<WasmModuleAny>;
    amrDetector: WasmModuleAny | null;
    wasmMemory: WebAssembly.Memory | null = null;
    amrWasmMemory: WebAssembly.Memory | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
        this.OrphosData = null;
        this.amrDetector = null;
        this.wasm = null;
        this.amrWasm = null;
        this.amrDetector = null;
        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg_orphos-bridge")
                .then((w) => {
                    this.wasm = w;
                    resolve(w);
                });
        });
        this.amrWasmPromise = new Promise((resolve) => {
            import("@/pkg_amr")
                .then((w) => {
                    this.amrWasm = w;
                    if (this.amrWasm.init_panic_hook) {
                        this.amrWasm.init_panic_hook();
                    }
                    resolve(w);
                });
        });
        import("@/pkg_orphos-bridge/index_bg.wasm").then((m) => { this.wasmMemory = m.memory; });
        import("@/pkg_amr/index_bg.wasm").then((m) => { this.amrWasmMemory = m.memory; });
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        if (!this.wasmMemory) return undefined;
        return this.wasmMemory.buffer.byteLength + (this.amrWasmMemory?.buffer.byteLength ?? 0);
    }

    step(s: string, fileName: string): void {
        this.worker.postMessage({ geneCallingStep: s, fileName });
    }

    waitForAmrWasm(): Promise<WasmModuleAny> {
        return this.amrWasm ? Promise.resolve(this.amrWasm) : this.amrWasmPromise;
    }

    async ensureAmrDetector(): Promise<WasmModuleAny> {
        if (this.amrDetector !== null) return this.amrDetector;
        const wasm = await this.waitForAmrWasm();
        this.amrDetector = await createAmrDetector(wasm);
        return this.amrDetector;
    }

    async callGenes(fileName: string, input_file: File, metag: boolean, closed_ends: boolean, mask: boolean, tt: number, non_sd: boolean, min_gene_fraction: number, min_gene_group_fraction: number): Promise<void> {
        console.log("Starting gene calling for: " + fileName);
        await this.waitForWasm();

        try {
            this.step('creating_interface', fileName);
            this.OrphosData = this.wasm.OrphosData.new(metag, "gff", closed_ends, mask, non_sd, tt);

            this.step('reading_fasta', fileName);
            this.OrphosData!.read_fasta(input_file);

            this.step('compressing_indexing', fileName);
            this.OrphosData!.index_fasta();

            this.step('calling_genes', fileName);
            this.OrphosData!.call_genes();
            this.step('genes_called', fileName);

            const geneMetadata = JSON.parse(this.OrphosData!.get_gene_metadata_json()) as GeneMetadataMap;
            let amrResult: AmrDetectionResult | null = null;
            let amrTsv = "";
            let amrError: string | null = null;
            let results: CallResult;
            try {
                this.step('detecting_amr', fileName);
                const cdsFasta = this.OrphosData!.get_cds_fasta();
                if (cdsFasta.trim().length > 0) {
                    const detector = await this.ensureAmrDetector();
                    const cdsBytes = new TextEncoder().encode(cdsFasta);
                    const amrJson = detector.detect_cds(
                        fileName,
                        cdsBytes,
                        min_gene_fraction,
                        min_gene_group_fraction
                    );
                    amrResult = JSON.parse(amrJson) as AmrDetectionResult;
                    amrTsv = buildGeneCallingAmrTsv(amrResult, geneMetadata);
                    this.step('writing_annotated_gff', fileName);
                    results = JSON.parse(this.OrphosData!.get_annotated_results("gff", amrJson));
                } else {
                    this.step('writing_gff', fileName);
                    results = JSON.parse(this.OrphosData!.get_results("gff"));
                }
            } catch (error) {
                amrError = error instanceof Error ? error.message : String(error);
                this.step('writing_gff', fileName);
                results = JSON.parse(this.OrphosData!.get_results("gff"));
            }
            this.step('done', fileName);

            const fastaBgz = this.OrphosData!.take_fasta_bgz();
            const fastaFai = this.OrphosData!.take_fasta_fai();
            const fastaGzi = this.OrphosData!.take_fasta_gzi();
            const gffBgz   = this.OrphosData!.take_gff_bgz();
            const gffCsi   = this.OrphosData!.take_gff_csi();

            this.worker.postMessage({
                fileName,
                output_file:    results.output_file,
                gene_count:     results.gene_count,
                sequence_count: results.sequence_count,
                fasta_bgz: fastaBgz,
                fasta_fai: fastaFai,
                fasta_gzi: fastaGzi,
                gff_bgz:   gffBgz,
                gff_csi:   gffCsi,
                gene_metadata: geneMetadata,
                amr_result: amrResult,
                amr_hit_count: amrResult ? amrResult.hits.length : 0,
                amr_tsv: amrTsv,
                amr_error: amrError,
                wasm_memory_bytes: this.memoryBytes(),
            }, [fastaBgz.buffer, fastaFai.buffer, fastaGzi.buffer, gffBgz.buffer, gffCsi.buffer]);
        } catch (error) {
            const message = error instanceof Error ? error.message : String(error);
            this.worker.postMessage({ error: true, fileName, message: message || 'unknown' });
        }
    }

    resetAll(): void {
        this.OrphosData = null;
        this.amrDetector = null;
    }
}
