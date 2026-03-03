import { Assembler } from './Assembler';

interface PreprocessMessage {
    preprocess: boolean;
    file1: File;
    file2: File | null;
    k: number;
    verbose: boolean;
    min_count: number;
    min_qual: number;
    csize: number;
    do_bloom: boolean;
    do_fit: boolean;
    no_bubble_collapse: boolean;
    no_dead_end_removal: boolean;
    use_gpu: boolean;
    gpu_power_pref: number;
}

interface AssembleMessage {
    assemble: boolean;
}

interface ResetMessage {
    reset: boolean;
}

interface ListAdaptersMessage {
    listAdapters: boolean;
}

type WorkerMessage = PreprocessMessage | AssembleMessage | ResetMessage | ListAdaptersMessage;

const ctx: Worker = self as unknown as Worker;
const assembler = new Assembler(ctx);

ctx.onmessage = async (evt: MessageEvent<WorkerMessage>) => {
    if (evt.data instanceof Object) {
        if ('preprocess' in evt.data && evt.data.preprocess) {
            console.log("Trying to preprocess!");
            const data = evt.data as PreprocessMessage;
            await assembler.preprocess(
                data.file1,
                data.file2,
                data.k,
                data.verbose,
                data.min_count,
                data.min_qual,
                data.csize,
                data.do_bloom,
                data.do_fit,
                data.no_bubble_collapse,
                data.no_dead_end_removal,
                data.use_gpu,
                data.gpu_power_pref,
            );
        } else if ('assemble' in evt.data && evt.data.assemble) {
            console.log("Trying to assemble!");
            await assembler.assemble();
        } else if ('reset' in evt.data && evt.data.reset) {
            console.log("Trying to reset!");
            assembler.resetAll();
        } else if ('listAdapters' in evt.data && evt.data.listAdapters) {
            console.log("Listing GPU adapters!");
            try {
                const adapters = await assembler.listAdapters();
                ctx.postMessage({ adaptersResult: adapters });
            } catch (error) {
                console.error("[Sparrowhawk] Adapter listing failed in worker:", error);
                ctx.postMessage({ adaptersResult: [] });
            }
        } else {
            throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
        }
    }
};
