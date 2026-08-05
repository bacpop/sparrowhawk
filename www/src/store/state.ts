import {AllResults, AllResultsSka, AllResultsSketchlib, AllResultsOrphos, AllResultsDeacon, AllResultsAmr, AllResultsEsm, EsmRetry, GpuAdapterInfo, TransmissionStandaloneResults, ProcessingState, ReadsPreprocessing, WorkerState,} from "@/types";


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
    allResults_amr: AllResultsAmr,
    allResults_esm: AllResultsEsm,
    esmRetry: EsmRetry,
    gpuAdapters: GpuAdapterInfo[],
    transmissionStandalone: TransmissionStandaloneResults,
    processingState: ProcessingState,
}
