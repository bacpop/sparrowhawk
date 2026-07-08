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

interface AlignMessage {
    align: boolean;
    files: File[];
    proportion_reads: number;
    k: number;
    rc: boolean;
    min_count: number;
    min_qual: number;
    qual_filter: number;
}

interface ClusterMessage {
    cluster: boolean;
    snp_threshold: number;
}

interface TransmissionClusterMessage {
    transmission_cluster: boolean;
    file: File;
    snp_threshold: number;
}

interface ResetMessage {
    reset: boolean;
}

type WorkerMessage = RefMessage | MapMessage | AlignMessage | ClusterMessage | TransmissionClusterMessage | ResetMessage;

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
        } else if ('align' in evt.data && evt.data.align) {
            const data = evt.data as AlignMessage;
            mapper.align(data.files, data.proportion_reads, data.rc, data.k, data.min_count, data.min_qual, data.qual_filter);
        } else if ('cluster' in evt.data && evt.data.cluster) {
            const data = evt.data as ClusterMessage;
            mapper.cluster(data.snp_threshold);
        } else if ('transmission_cluster' in evt.data && evt.data.transmission_cluster) {
            const data = evt.data as TransmissionClusterMessage;
            mapper.clusterUploadedAlignment(data.file, data.snp_threshold);
        } else if ('reset' in evt.data && evt.data.reset) {
            mapper.resetAll();
        } else {
            throw new Error("Event " + JSON.stringify(evt.data) + " is not supported");
        }
    }
};
