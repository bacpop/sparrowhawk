import { Embedder } from './Embedder';

interface EmbedMessage {
    embed: true;
    file: File;
    sampleName: string;
    use_gpu: boolean;
    gpu_power_pref: number;
    gpu_tasks_max: number;
}

interface ResetMessage {
    reset: true;
}

type WorkerMessage = EmbedMessage | ResetMessage;

const ctx: Worker = self as unknown as Worker;
const embedder = new Embedder(ctx);

// For catching errors from Rust, forcing also to get again a new worker.
function reportFatal(detail: unknown): void {
    const message = detail instanceof Error ? detail.message : String(detail);
    console.error("[Sparrowhawk] embedder aborted:", message);
    ctx.postMessage({
        error: true,
        wasmPanic: true,
        wasGpu: embedder.isUsingGpu(),
        sampleName: embedder.getCurrentSample(),
        message,
    });
}

// These "addEventListener" allow to catch errors this way, this might be useful for the other workers in the future
self.addEventListener("error", (event) => {
    reportFatal((event as ErrorEvent).error ?? (event as ErrorEvent).message);
});
self.addEventListener("unhandledrejection", (event) => {
    reportFatal((event as PromiseRejectionEvent).reason);
});


ctx.onmessage = async (evt: MessageEvent<WorkerMessage>) => {
    if (!(evt.data instanceof Object)) return;

    if ('embed' in evt.data && evt.data.embed) {
        const data = evt.data as EmbedMessage;
        await embedder.embedThisFile(data.file, data.sampleName, {
            use_gpu: data.use_gpu,
            gpu_power_pref: data.gpu_power_pref,
            gpu_tasks_max: data.gpu_tasks_max ?? 0,
        });
    } else if ('reset' in evt.data && evt.data.reset) {
        embedder.resetAll();
    } else {
        throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
    }
};
