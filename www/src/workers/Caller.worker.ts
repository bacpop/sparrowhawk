import { Caller } from './Caller';

interface CallMessage {
    call: boolean;
    fileName: string;
    input_file: File;
    metag: boolean;
    closed_ends: boolean;
    mask: boolean;
    tt: number;
    non_sd: boolean;
    min_gene_fraction: number;
    min_gene_group_fraction: number;
}

interface ResetMessage {
    reset: boolean;
}

type WorkerMessage = CallMessage | ResetMessage;

const ctx: Worker = self as unknown as Worker;
const caller = new Caller(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (evt.data instanceof Object) {
        if ('call' in evt.data && evt.data.call) {
            const data = evt.data as CallMessage;
            caller.callGenes(data.fileName, data.input_file, data.metag, data.closed_ends, data.mask, data.tt, data.non_sd, data.min_gene_fraction, data.min_gene_group_fraction);
        } else if ('reset' in evt.data && evt.data.reset) {
            caller.resetAll();
        } else {
            throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
        }
    }
};