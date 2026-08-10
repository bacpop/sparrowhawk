interface MapResult {
    "Number of variants": number;
    "Coverage": number;
    "Mapped sequences": string[];
    "VCF": string;
}

interface AlignResult {
    names: string[];
    newick: string;
    alignment: string;
}

interface ClusterLabels {
    [sampleName: string]: number;
}

interface TransmissionGraphNode { id: string; cluster: number; }
interface TransmissionGraphLink { source: string; target: string; snp_distance: number; }
interface TransmissionGraphData {
    nodes: TransmissionGraphNode[];
    links: TransmissionGraphLink[];
}

interface SkaData {
    get_reference(): string;
    map(file: File, revReadFile: File | null, proportion_reads: number, min_count: number, min_qual: number, qual_filter: number, sampleName: string): string;
}

interface AlignData {
    align(files: File[], proportion_reads: number, min_count: number, min_qual: number, qual_filter: number): string;
    get_graph_json(snp_threshold: number): string;
    get_distances_csv(): string;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

export class Mapper {
    worker: Worker;
    wasm: WasmModuleAny | null;
    SkaData: SkaData | null;
    AlignData: AlignData | null;
    wasmPromise: Promise<WasmModuleAny>;
    wasmMemory: WebAssembly.Memory | null = null;

    constructor(worker: Worker) {
        this.worker = worker;
        this.SkaData = null;
        this.AlignData = null;
        this.wasm = null;
        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg_ska")
                .then((w) => {
                    this.wasm = w;
                    resolve(w);
                });
        });
        import("@/pkg_ska/index_bg.wasm").then((m) => { this.wasmMemory = m.memory; });
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    memoryBytes(): number | undefined {
        return this.wasmMemory ? this.wasmMemory.buffer.byteLength : undefined;
    }

    async set_ref(file: File, k: number, rc: boolean, ambig_mask: boolean, repeat_mask: boolean): Promise<void> {
        await this.waitForWasm();

        if (this.SkaData === null) {
            this.SkaData = this.wasm.SkaData.new(file, k, rc, ambig_mask, repeat_mask);
        }
        this.worker.postMessage({
            ref: file,
            sequences: this.SkaData!.get_reference().split('\n')
        });
    }

    map(file: File, revReadFile: File | null, proportion_reads: number, min_count: number, min_qual: number, qual_filter: number): void {
        console.log("Mapping reads to reference with proportion_reads: " + proportion_reads);
        if (this.SkaData === null) {
            throw new Error("SkaRef::map - reference does not exist yet.");
        }

        try {
            const t0 = performance.now();
            const outname = (revReadFile != null) ? file.name.replace(new RegExp("(?:_1)?\\.(?:fa|fna|fasta|fq|fnq|fastq)(?:\\.gz)?" + String.fromCharCode(36)), "") : file.name;
            const results: MapResult = JSON.parse(this.SkaData.map(file, revReadFile, proportion_reads, min_count, min_qual, qual_filter, outname));

            this.worker.postMessage({
                nb_variants: results["Number of variants"],
                coverage: results["Coverage"],
                name: outname,
                mapped_sequences: results["Mapped sequences"],
                mapping_vcf: results["VCF"],
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes(),
            });
        } catch {
            this.worker.postMessage({ error: true, message: 'memory' });
        }
    }

    align(files: File[], proportion_reads: number, rc: boolean, k: number, min_count: number, min_qual: number, qual_filter: number): void {
        console.log("Processing uploaded fastX files with proportion_reads: " + proportion_reads + " and k: " + k);

        try {
            const t0 = performance.now();
            if (this.AlignData === null) {
                this.AlignData = this.wasm.AlignData.new(k, rc);
            }

            const results: AlignResult = JSON.parse(this.AlignData!.align(files, proportion_reads, min_count, min_qual, qual_filter));

            const distances_csv = this.AlignData!.get_distances_csv();
            this.worker.postMessage({
                aligned: true,
                names: results.names,
                newick: results.newick,
                alignment: results.alignment,
                distances_csv,
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes(),
            });
        } catch {
            this.worker.postMessage({ error: true, message: 'memory' });
        }
    }

    async clusterUploadedAlignment(file: File, snp_threshold: number): Promise<void> {
        console.log("Running standalone transmission clustering with SNP threshold: " + snp_threshold);

        try {
            const t0 = performance.now();
            await this.waitForWasm();
            const alignmentText = await file.text();
            const alignData = this.wasm.AlignData.from_alignment_text(alignmentText);
            const clusters: ClusterLabels = JSON.parse(this.wasm.ska_cluster(alignData, snp_threshold));
            const graph: TransmissionGraphData = JSON.parse(alignData.get_graph_json(snp_threshold));
            this.worker.postMessage({ clustered: true, clusters, graph,
                elapsedMs: Math.round(performance.now() - t0),
                wasmMemoryBytes: this.memoryBytes() });
        } catch (error) {
            const message = error instanceof Error ? error.message : String(error);
            this.worker.postMessage({ error: true, message });
        }
    }

    cluster(snp_threshold: number): void {
        console.log("Running transmission clustering with SNP threshold: " + snp_threshold);
        if (this.AlignData === null) {
            throw new Error("Mapper::cluster - no alignment data available.");
        }

        try {
            const clusters: ClusterLabels = JSON.parse(this.wasm.ska_cluster(this.AlignData, snp_threshold));
            const graph: TransmissionGraphData = JSON.parse(this.AlignData.get_graph_json(snp_threshold));
            this.worker.postMessage({ clustered: true, clusters, graph });
        } catch {
            this.worker.postMessage({ error: true, message: 'memory' });
        }
    }

    resetAll(): void {
        this.SkaData = null;
        this.AlignData = null;
    }
}
