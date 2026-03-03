import {AllResults, AllResultsSka, AllResultsSketchlib, AllResultsOrphos, AllResultsDeacon, ProcessingState, ReadsPreprocessing, WorkerState,} from "@/types";


export interface RootState {
    workerState: WorkerState,
    readsFileNames: string | null,
    readsPreprocessing: ReadsPreprocessing,
    allResults: AllResults,
    refSet: string | null,
    allResults_ska: AllResultsSka,
    errors: string,
    min_count: number,
    allResults_sketchlib: AllResultsSketchlib,
    allResults_orphos: AllResultsOrphos,
    allResults_deacon: AllResultsDeacon,
    processingState: ProcessingState,
    gpuAdapters: { index: number; name: string }[],
}
