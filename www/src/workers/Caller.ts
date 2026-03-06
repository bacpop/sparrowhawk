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

    constructor(worker: Worker) {
        this.worker = worker;
        this.OrphosData = null;
        this.wasm = null;
        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg_orphos-bridge")
                .then((w) => {
                    this.wasm = w;
                    resolve(w);
                });
        });
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    step(s: string, fileName: string): void {
        this.worker.postMessage({ geneCallingStep: s, fileName });
    }

    async callGenes(fileName: string, input_file: File, metag: boolean, closed_ends: boolean, mask: boolean, tt: number, non_sd: boolean): Promise<void> {
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

            this.step('writing_gff', fileName);
            const results: CallResult = JSON.parse(this.OrphosData!.get_results("gff"));
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
            }, [fastaBgz.buffer, fastaFai.buffer, fastaGzi.buffer, gffBgz.buffer, gffCsi.buffer]);
        } catch {
            this.worker.postMessage({ error: true, fileName, message: 'memory' });
        }
    }

    resetAll(): void {
        this.OrphosData = null;
    }
}
