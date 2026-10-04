import { Mapper } from './Mapper';

interface RefMessage {
    ref: boolean;
    file: File;
    k: number;
    rc: boolean;
    ambig_mask: boolean;
    repeat_mask: boolean;
}

interface MapMessage {
    map: boolean;
    file: File;
    revReads: File | null;
    proportion_reads: number;
    min_count: number;
    min_qual: number;
    qual_filter: number;
}

interface AlignmentStartMessage {
    append: boolean;
    rc: boolean;
    alignStart: boolean;
    runId: number;
    fileNames: string[];
    k: number;
}

interface AlignmentSampleMessage {
    addAlignmentSample: boolean;
    runId: number;
    sampleIndex: number;
    sampleName: string;
    packedKmers: Uint32Array;
}

interface AlignmentFinishMessage {
    finishAlignment: boolean;
    runId: number;
}

interface AlignmentCancelMessage {
    cancelAlignment: boolean;
    runId: number;
}

interface AlignmentExportMessage {
    exportAlignment: boolean;
    runId: number;
    requestId: number;
}

interface ClusterMessage {
    runId: number;
    requestId: number;
    cluster: boolean;
    snp_threshold: number;
}

interface TransmissionClusterMessage {
    requestId: number;
    transmission_cluster: boolean;
    file: File;
    snp_threshold: number;
}

interface ResetMessage {
    reset: boolean;
}

type WorkerMessage = RefMessage | MapMessage | AlignmentStartMessage | AlignmentSampleMessage |
    AlignmentFinishMessage | AlignmentCancelMessage | ClusterMessage | TransmissionClusterMessage |
    AlignmentExportMessage | ResetMessage;

const ctx: Worker = self as unknown as Worker;
const mapper = new Mapper(ctx);

ctx.onmessage = (evt: MessageEvent<WorkerMessage>) => {
    if (evt.data instanceof Object) {
        if ('ref' in evt.data && evt.data.ref) {
            const data = evt.data as RefMessage;
            mapper.set_ref(data.file, data.k, data.rc, data.ambig_mask, data.repeat_mask);
        } else if ('map' in evt.data && evt.data.map) {
            const data = evt.data as MapMessage;
            mapper.map(data.file, data.revReads, data.proportion_reads, data.min_count, data.min_qual, data.qual_filter);
        } else if ('alignStart' in evt.data && evt.data.alignStart) {
            const data = evt.data as AlignmentStartMessage;
            void mapper.beginAlignment(data.runId, data.fileNames, data.k, data.rc, data.append);
        } else if ('addAlignmentSample' in evt.data && evt.data.addAlignmentSample) {
            const data = evt.data as AlignmentSampleMessage;
            mapper.addAlignmentSample(data.runId, data.sampleIndex, data.sampleName, data.packedKmers);
        } else if ('finishAlignment' in evt.data && evt.data.finishAlignment) {
            const data = evt.data as AlignmentFinishMessage;
            mapper.finishAlignment(data.runId);
        } else if ('cancelAlignment' in evt.data && evt.data.cancelAlignment) {
            const data = evt.data as AlignmentCancelMessage;
            mapper.cancelAlignment(data.runId);
        } else if ('exportAlignment' in evt.data && evt.data.exportAlignment) {
            const data = evt.data as AlignmentExportMessage;
            mapper.exportAlignment(data.runId, data.requestId);
        } else if ('cluster' in evt.data && evt.data.cluster) {
            const data = evt.data as ClusterMessage;
            mapper.cluster(data.snp_threshold, data.runId, data.requestId);
        } else if ('transmission_cluster' in evt.data && evt.data.transmission_cluster) {
            const data = evt.data as TransmissionClusterMessage;
            void mapper.clusterUploadedAlignment(data.file, data.snp_threshold, data.requestId);
        } else if ('reset' in evt.data && evt.data.reset) {
            mapper.resetAll();
        } else {
            throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
        }
    }
};
