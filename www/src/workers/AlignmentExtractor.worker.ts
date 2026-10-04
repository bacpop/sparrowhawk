import { loadWasmPackage } from "@/wasm/loader";
import { skaPackage } from "@/wasm/packages/ska";

interface ExtractSampleMessage {
    extractSample: true;
    runId: number;
    sampleIndex: number;
    sampleName: string;
    firstFile: File;
    secondFile: File | null;
    k: number;
    rc: boolean;
    proportion_reads: number;
    min_count: number;
    min_qual: number;
    qual_filter: number;
}

type SkaModule = typeof import("@/pkg_ska");

const ctx: Worker = self as unknown as Worker;
let wasmPromise: Promise<SkaModule> | null = null;

function loadWasm(): Promise<SkaModule> {
    if (wasmPromise === null) {
        wasmPromise = loadWasmPackage(skaPackage, (notice) => ctx.postMessage(notice))
            .then((loaded) => loaded.api);
        void wasmPromise.catch(() => { wasmPromise = null; });
    }
    return wasmPromise;
}

ctx.onmessage = async (evt: MessageEvent<ExtractSampleMessage>) => {
    const data = evt.data;
    if (!(data instanceof Object) || !data.extractSample) return;

    const startedAt = performance.now();
    try {
        const wasm = await loadWasm();
        const packedKmers: Uint32Array = wasm.extract_alignment_sample(
            data.firstFile,
            data.secondFile,
            data.k,
            data.sampleIndex,
            data.sampleName,
            data.rc,
            data.proportion_reads,
            data.min_count,
            data.min_qual,
            data.qual_filter,
        );
        ctx.postMessage({
            type: "sampleExtracted",
            runId: data.runId,
            sampleIndex: data.sampleIndex,
            sampleName: data.sampleName,
            packedKmers,
            elapsedMs: Math.round(performance.now() - startedAt),
        }, [packedKmers.buffer as ArrayBuffer]);
    } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Sample extraction failed", {
            runId: data.runId,
            sampleIndex: data.sampleIndex,
            sampleName: data.sampleName,
            error,
        });
        ctx.postMessage({
            type: "sampleError",
            runId: data.runId,
            sampleIndex: data.sampleIndex,
            sampleName: data.sampleName,
            message,
        });
    }
};
