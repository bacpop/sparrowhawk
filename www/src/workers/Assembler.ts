interface PreprocessingResult {
    nkmers: number;
    histo: number[];
    used_min_count: number;
}

interface AssemblyResult {
    outfasta: string;
    ncontigs: number;
    outdot: string;
    outgfa: string;
    outgfav2: string;
}

interface WasmModule {
    AssemblyHelper: {
        new(
            k: number,
            verbose: boolean,
            min_count: number,
            min_qual: number,
            chunk_size: number,
            do_bloom: boolean,
            do_fit: boolean,
            no_bubble_collapse: boolean,
            no_dead_end_removal: boolean,
            use_gpu: boolean,
            gpu_power_pref: number,
        ): AssemblyHelper;
    };
    list_gpu_adapters(): Promise<string>; // returns JSON string
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type WasmModuleAny = any;

interface AssemblyHelper {
    preprocess(file1: File, file2: File | null): Promise<void>;
    get_preprocessing_info(): string;
    assemble(): void;
    get_assembly(): string;
}

function getWebGLRendererLabel(): string {
    try {
        const canvas = new OffscreenCanvas(1, 1);
        const gl: any = canvas.getContext('webgl2') || canvas.getContext('webgl');
        if (!gl) return '';
        const ext = gl.getExtension('WEBGL_debug_renderer_info');
        if (!ext) return '';
        let renderer: string = gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) || '';
        if (!renderer) return '';

        // Extract the GPU model from the ANGLE wrapper when present:
        //   "ANGLE (Vendor, <model>, Backend)" → "<model>"
        const m = renderer.match(/^ANGLE\s*\([^,]+,\s*(.+),\s*[^,]+\)$/);
        if (m) renderer = m[1];

        // Strip backend-description suffixes that describe rendering API, not the GPU:
        //   macOS:   "ANGLE Metal Renderer: Apple M1 Pro"  → "Apple M1 Pro"
        //   Windows: "Intel(R) Iris(R) Xe Graphics Direct3D11 vs_5_0 ps_5_0" → "Intel(R) Iris(R) Xe Graphics"
        renderer = renderer
            .replace(/^ANGLE\s+Metal\s+Renderer:\s*/i, '')
            .replace(/\s+Direct3D\d+(\s+vs_\d+_\d+\s+ps_\d+_\d+)?$/i, '')
            .trim();

        return renderer;
    } catch {
        return '';
    }
}

export class Assembler {
    worker: Worker;
    wasm: WasmModuleAny | null;
    helper: AssemblyHelper | null;
    wasmPromise: Promise<WasmModuleAny>;
    noBubbleCollapse: boolean = false;
    noDeadEndRemoval: boolean = false;

    constructor(worker: Worker) {
        this.worker = worker;
        this.helper = null;
        this.wasm = null;

        this.wasmPromise = new Promise((resolve) => {
            import("@/pkg")
                .then((w) => {
                    this.wasm = w;
                    resolve(w);
                });
        });
    }

    waitForWasm(): Promise<WasmModuleAny> {
        return this.wasm ? Promise.resolve(this.wasm) : this.wasmPromise;
    }

    async preprocess(
        file1: File,
        file2: File | null,
        k: number,
        verbose: boolean,
        min_count: number,
        min_qual: number,
        csize: number,
        do_bloom: boolean,
        do_fit: boolean,
        no_bubble_collapse: boolean,
        no_dead_end_removal: boolean,
        use_gpu: boolean,
        gpu_power_pref: number,
    ): Promise<void> {
        await this.waitForWasm();

        this.noBubbleCollapse = no_bubble_collapse;
        this.noDeadEndRemoval = no_dead_end_removal;

        if (this.helper === null) {
            try {
                this.helper = this.wasm!.AssemblyHelper.new(
                    k, verbose,
                    min_count, min_qual,
                    csize, do_bloom, do_fit,
                    no_bubble_collapse, no_dead_end_removal,
                    use_gpu, gpu_power_pref,
                );
                this.helper!.preprocess(file1, file2);
            } catch (error) {
                console.log("Webassembly error found! Most surely, memory issue.");
                console.error(error);
                this.worker.postMessage({ reset: true });
                return;
            }
        }

        await this.helper!.preprocess(file1, file2);

        const resultsjson: PreprocessingResult = JSON.parse(
            this.helper!.get_preprocessing_info()
        );

        this.worker.postMessage({
            nKmers:         resultsjson.nkmers,
            histo:          resultsjson.histo,
            used_min_count: resultsjson.used_min_count,
        });
    }

    async assemble(): Promise<void> {
        console.log("Initiating assembly from worker");
        this.helper!.assemble();

        console.log("Assembly finished");
        const resultsjson: AssemblyResult = JSON.parse(this.helper!.get_assembly());

        console.log("Posting output as message");
        this.worker.postMessage({
            outfasta: resultsjson.outfasta,
            ncontigs: resultsjson.ncontigs,
            outdot: resultsjson.outdot,
            outgfa: resultsjson.outgfa,
            outgfav2: resultsjson.outgfav2
        });
    }

    async listAdapters(): Promise<{index: number, name: string}[]> {
        // Use the native WebGPU JS API directly — wgpu's request_adapter() can hang
        // in some browser/driver combinations. The JS API is reliable and doesn't
        // require a wasm rebuild. Indices match the Rust gpu_power_pref mapping:
        //   0 → None (browser decides), 1 → HighPerformance, 2 → LowPower
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const gpu = (navigator as any).gpu;
        if (!gpu) {
            console.warn("[Sparrowhawk] WebGPU not available (navigator.gpu is undefined)");
            return [];
        }

        const webglRenderer = getWebGLRendererLabel();

        // Pass 1: collect (index, hwLabel) for every power-preference that returns an adapter.
        interface AdapterEntry { index: number; hwLabel: string; }
        const entries: AdapterEntry[] = [];

        const prefs: [number, string | undefined][] = [
            [0, undefined],
            [1, 'high-performance'],
            [2, 'low-power'],
        ];

        for (const [idx, powerPreference] of prefs) {
            try {
                // eslint-disable-next-line @typescript-eslint/no-explicit-any
                const adapter: any = await gpu.requestAdapter(
                    powerPreference !== undefined ? { powerPreference } : undefined
                );
                if (adapter) {
                    const info = adapter.info;
                    const description = (info?.description || '').trim();
                    const vendorArch = (info?.vendor && info?.architecture)
                        ? `${info.vendor} (${info.architecture})`
                        : '';
                    const vendorDevice = info?.vendor
                        ? `${info.vendor}${info?.device ? ` ${info.device}` : ''}`
                        : (info?.device || '');
                    const hwLabel = description || vendorArch || vendorDevice || webglRenderer || 'GPU';
                    entries.push({ index: idx, hwLabel });
                }
            } catch (e) {
                console.warn(`[Sparrowhawk] requestAdapter(${powerPreference ?? 'default'}) failed:`, e);
            }
        }

        // Pass 2: build result list.
        // If only one distinct GPU exists → one clean entry (no mode suffix), prefer index 1.
        // If multiple distinct GPUs → include mode suffix so the user can disambiguate.
        const uniqueHwLabels = new Set(entries.map(e => e.hwLabel));
        const modeLabels: Record<number, string> = {
            0: 'default',
            1: 'high-performance',
            2: 'low-power',
        };

        const results: {index: number, name: string}[] = [];
        const seen = new Set<string>();

        if (uniqueHwLabels.size <= 1) {
            // Single GPU — show once with hardware name only.
            // Prefer index 1 (high-performance); fall back to whatever index was found.
            const preferred = entries.find(e => e.index === 1) ?? entries[0];
            if (preferred) {
                results.push({ index: preferred.index, name: preferred.hwLabel });
            }
        } else {
            // Multiple distinct GPUs — suffix with mode for disambiguation.
            for (const { index, hwLabel } of entries) {
                const modeLabel = modeLabels[index] ?? String(index);
                const name = `${hwLabel} — ${modeLabel}`;
                if (!seen.has(name)) {
                    seen.add(name);
                    results.push({ index, name });
                }
            }
        }

        console.log("[Sparrowhawk] Detected GPU adapters:", results);
        return results;
    }

    resetAll(): void {
        this.helper = null;
    }
}
