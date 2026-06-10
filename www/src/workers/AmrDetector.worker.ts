import { AmrDetectorWorker } from './AmrDetector';

interface DetectAmrMessage {
    detectAmr: true;
    file: File;
    sampleName: string;
    min_gene_fraction: number;
    min_gene_group_fraction: number;
}

type WorkerMessage = DetectAmrMessage;

const ctx: Worker = self as unknown as Worker;
const detector = new AmrDetectorWorker(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (!(evt.data instanceof Object)) return;
    if ('detectAmr' in evt.data && evt.data.detectAmr) {
        const data = evt.data as DetectAmrMessage;
        detector.detectThisFile(
            data.file,
            data.sampleName,
            data.min_gene_fraction,
            data.min_gene_group_fraction
        );
    } else {
        throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
    }
};
