import { Assembler } from './Assembler';

interface PreprocessMessage {
    preprocess: boolean;
    file1: File;
    file2: File | null;
    k: number;
    verbose: boolean;
    min_count: number;
    min_qual: number | null;
    do_bloom: boolean;
    do_fit: boolean;
    multik: boolean;
    ks: number[];
    min_contig_length: number;
}

interface AssembleMessage {
    assemble: boolean;
}

interface ResetMessage {
    reset: boolean;
}

type WorkerMessage = PreprocessMessage | AssembleMessage | ResetMessage;

const ctx: Worker = self as unknown as Worker;
const assembler = new Assembler(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (evt.data instanceof Object) {
        if ('preprocess' in evt.data && evt.data.preprocess) {
            console.log("Trying to preprocess!");
            const data = evt.data as PreprocessMessage;
            assembler.preprocess(
                data.file1,
                data.file2,
                data.k,
                data.verbose,
                data.min_count,
                data.min_qual,
                data.do_bloom,
                data.do_fit,
                data.multik,
                data.ks,
                data.min_contig_length
            );
        } else if ('assemble' in evt.data && evt.data.assemble) {
            console.log("Trying to assemble!");
            assembler.assemble();
        } else if ('reset' in evt.data && evt.data.reset) {
            console.log("Trying to reset!");
            assembler.resetAll();
        } else {
            throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
        }
    }
};
