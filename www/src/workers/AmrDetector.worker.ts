import { AmrDetector } from "./AmrDetector";

interface LoadIndexMessage { loadIndex: true; file: File; }
interface DetectMessage { detect: true; file: File; sampleName: string; min_gene_hits: number; min_family_hits: number; }
interface ResetMessage { reset: true; }
type WorkerMessage = LoadIndexMessage | DetectMessage | ResetMessage;

const ctx: Worker = self as unknown as Worker;
const detector = new AmrDetector(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (!(evt.data instanceof Object)) return;
    if ("loadIndex" in evt.data && evt.data.loadIndex) {
        detector.loadIndex(evt.data.file);
    } else if ("detect" in evt.data && evt.data.detect) {
        detector.detect(evt.data.file, evt.data.sampleName, evt.data.min_gene_hits, evt.data.min_family_hits);
    } else if ("reset" in evt.data && evt.data.reset) {
        detector.resetAll();
    }
};
