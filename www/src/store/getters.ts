import {RootState} from "@/store/state";
import {GetterTree} from "vuex";

export const getters: GetterTree<RootState, RootState> = {
    // Processing state getters
    isPreprocessing(state: RootState): boolean {
        return state.processingState.isPreprocessing;
    },
    isAssembling(state: RootState): boolean {
        return state.processingState.isAssembling;
    },
    isIndexingRef(state: RootState): boolean {
        return state.processingState.isIndexingRef;
    },
    isMapping(state: RootState): boolean {
        return state.processingState.isMapping;
    },
    isMappingFiles(state: RootState): Set<string> {
        return state.processingState.isMappingFiles;
    },
    isAligning(state: RootState): boolean {
        return state.processingState.isAligning;
    },
    isIdentifying(state: RootState): boolean {
        return state.processingState.isIdentifying;
    },
    isAnyProcessing(state: RootState): boolean {
        return state.processingState.isPreprocessing ||
            state.processingState.isAssembling ||
            state.processingState.isIndexingRef ||
            state.processingState.isMapping ||
            state.processingState.isAligning ||
            state.processingState.isIdentifying ||
            state.processingState.isCallingGenes ||
            state.processingState.isFilteringDeacon ||
            state.processingState.isDetectingAmr;
    },
    assemblyState(state: RootState): string {
        return state.processingState.assemblyState;
    },
    readsProcessing(state: RootState) {
        return (state.readsFileNames != null && state.allResults.nContigs == null);
    },
    readsPreprocessing(state: RootState) {
        return (state.readsFileNames != null && state.allResults.nContigs == null && state.readsPreprocessing.nKmers == null);
    },
    readsPreprocessed(state: RootState) {
        return (state.readsPreprocessing.nKmers != null);
    },
    queryAssembled(state: RootState) {
        return state.allResults.fastaOutput != "";
    },
    getErrors(state: RootState) {
        return (state.errors != "");
    },
    getErrorType(state: RootState) {
        return state.errors;
    },
    readsName(state: RootState) {
        return state.readsFileNames;
    },

    // SKA
    refProcessed(state: RootState) {
        console.log("Ref processed (getters)")
        return state.refSet != null;
    },
    refName(state: RootState) {
        return state.refSet == null ? "" : state.refSet;
    },
    queryProcessed(state: RootState) {
        const first_input: string = Object.keys(state.allResults_ska.mapResults)[0]
        console.log(state.allResults_ska.mapResults[first_input] ? state.allResults_ska.mapResults[first_input].mapped : false)
        return state.allResults_ska.mapResults[first_input] ? state.allResults_ska.mapResults[first_input].mapped : false;
    },
    alignmentProcessed(state: RootState) {
        return state.allResults_ska.alignResults[0] ? state.allResults_ska.alignResults[0].aligned : false;
    },

    // SKETCHLIB
    sampleIdentified(state: RootState) {
        return Object.keys(state.allResults_sketchlib.results).length > 0;
    },
    isIdentifyingFiles(state: RootState): Set<string> {
        return state.processingState.isIdentifyingFiles;
    },
    sketchlibResults(state: RootState) {
        return state.allResults_sketchlib.results;
    },
    
    // ORPHOS
    genesCalled(state: RootState) {
        return Object.keys(state.allResults_orphos.results).length > 0;
    },
    callingGenes(state: RootState) {
        return state.processingState.isCallingGenes;
    },
    orphosResults(state: RootState) {
        return state.allResults_orphos.results;
    },
    callingGenesFiles(state: RootState): Set<string> {
        return state.processingState.isCallingGenesFiles;
    },
    geneCallingStep(state: RootState): string {
        return state.processingState.geneCallingStep;
    },
    geneCallingProgressTotal(state: RootState): number {
        return state.processingState.geneCallingProgressTotal;
    },

    // SKA clustering
    hasClusterResults(state: RootState): boolean { return state.allResults_ska.clusterResults !== null; },
    isClustering(state: RootState): boolean      { return state.processingState.isClustering; },

    // SKA errors
    skaError(state: RootState)       { return state.allResults_ska.error; },
    // SKETCHLIB errors
    sketchlibError(state: RootState) { return state.allResults_sketchlib.error; },
    // ORPHOS errors
    orphosError(state: RootState)    { return state.allResults_orphos.error; },
    // DEACON errors
    deaconError(state: RootState)    { return state.allResults_deacon.error; },
    amrError(state: RootState)       { return state.allResults_amr.error; },

    // DEACON
    deaconIndexLoaded(state: RootState)    { return state.allResults_deacon.indexLoaded; },
    isLoadingDeaconIndex(state: RootState) { return state.allResults_deacon.isLoadingIndex; },
    isFilteringDeacon(state: RootState)    { return state.processingState.isFilteringDeacon; },
    filteringDeaconFiles(state: RootState): Set<string> { return state.processingState.isFilteringDeaconFiles; },
    deaconFiltered(state: RootState)       { return Object.keys(state.allResults_deacon.results).length > 0; },
    deaconResults(state: RootState)        { return state.allResults_deacon.results; },
    amrIndexLoaded(state: RootState)       { return state.allResults_amr.indexLoaded; },
    isLoadingAmrIndex(state: RootState)    { return state.allResults_amr.isLoadingIndex; },
    isDetectingAmr(state: RootState)       { return state.processingState.isDetectingAmr; },
    detectingAmrFiles(state: RootState): Set<string> { return state.processingState.isDetectingAmrFiles; },
    amrDetected(state: RootState)          { return Object.keys(state.allResults_amr.results).length > 0; },
    amrResults(state: RootState)           { return state.allResults_amr.results; },
}
