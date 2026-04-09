import {RootState} from "@/store/state";

export const emptyState = (): RootState => ({
    readsFileNames: null,
    errors: "",
    min_count: 0,
    refSet: null,

    readsPreprocessing: {
        nKmers: null,
        histo: [],
        used_min_count: null,
    },

    allResults: {
        nContigs: null,
        fastaOutput: "",
        dotOutput: "",
        gfaOutput: "",
        gfav2Output: "",
    },

    allResults_ska: {
        mapResults: {},
        alignResults: {},
        ref: [],
        error: null,
        clusterResults: null,
        transmissionGraph: null,
    },

    allResults_sketchlib: {
        results: {},
        error: null,
    },

    allResults_orphos: {
        results: {},
        error: null,
    },

    allResults_deacon: {
        indexFileName: null,
        indexInfo: null,
        indexLoaded: false,
        isLoadingIndex: false,
        results: {},
        error: null,
    },

    workerState: {
        worker: null,
        worker_ska: null,
        worker_deacon: null,
        workers_sketchlib: [],
        workers_orphos: [],
    },

    processingState: {
        isPreprocessing: false,
        isAssembling: false,
        isIndexingRef: false,
        isMapping: false,
        isMappingFiles: new Set<string>(),
        isAligning: false,
        isIdentifying: false,
        isIdentifyingFiles: new Set<string>(),
        assemblyState: '',
        isCallingGenes: false,
        isCallingGenesFiles: new Set<string>(),
        geneCallingStep: '',
        geneCallingProgressTotal: 0,
        isFilteringDeacon: false,
        isFilteringDeaconFiles: new Set<string>(),
        isClustering: false,
    },
});

export const findReadPair = (fileName: string, files: Array<File>): {
    pairFile: File | undefined,
    sampleName: string
} => {
    const baseName = fileName.replace(regExpWithOneOnly, '');
    const pairFile = files.find(f => regExpWithTwoOnly.test(f.name));
    return {pairFile, sampleName: baseName};
};

export const fastaExtensionsList = ["fa", "fna", "fasta"];
export const fastqExtensionsList = ["fq", "fnq", "fastq"];

export const fastaExtensionsWithDotAndCompressList = fastaExtensionsList.map(x => "." + x).concat(fastaExtensionsList.map(x => "." + x + ".gz"));
export const fastqExtensionsWithDotAndCompressList = fastqExtensionsList.map(x => "." + x).concat(fastqExtensionsList.map(x => "." + x + ".gz"));
export const fastxExtensionsWithDotAndCompressList = fastaExtensionsWithDotAndCompressList.concat(fastqExtensionsWithDotAndCompressList);

export const fastaExtensionsReg = "\\.(?:" + fastaExtensionsList.join("|") + ")(?:\\.gz)?";
export const fastqExtensionsReg = "\\.(?:" + fastqExtensionsList.join("|") + ")(?:\\.gz)?";
export const fastxExtensionsReg = "\\.(?:" + fastaExtensionsList.concat(fastqExtensionsList).join("|") + ")(?:\\.gz)?";

export const regExpWithTwoNumbers = new RegExp("_(?:1|2)" + fastqExtensionsReg + "$");
export const regExpWithOneOnly = new RegExp("_1" + fastqExtensionsReg + "$");
export const regExpWithTwoOnly = new RegExp("_2" + fastqExtensionsReg + "$");
export const regExpForAnyFastx = new RegExp(fastxExtensionsReg + "$");
export const regExpForAnyFasta = new RegExp(fastaExtensionsReg + "$");

export const getFilesToProcess = (files: Array<File>) => {
    const outlist = [];
    const tmplist: string[] = [];
    // Global idea: test if it's fasta, if it is, then it is an individual file to process.
    // Then, if not, test if it ends in 1 or 2 before the extension. If it doesn't, then it's alone as well.
    // If it does, search for its pair. If it doesn't have one, it's alone. If it has one, it's together.
    // In this last case, make sure to add it to the tmplist when replacing its extension plus ID.
    console.log("Length of input list: " + files.length);
    // console.log(regExpForAnyFasta);
    // console.log(regExpWithTwoNumbers);
    for (let i = files.length - 1; i >= 0; i--) {
        if (regExpForAnyFasta.test(files[i].name)) {
            console.log("This file is FASTA");
            outlist.push([i]);
        } else {
            // Not FASTA
            // Does it end on 1 or 2 before extension?
            console.log("This file is not FASTA");
            if (regExpWithTwoNumbers.test(files[i].name)) {
                // Yes, now we have to look for its partner
                const tmpnam = files[i].name.replace(regExpWithTwoNumbers, "");
                // check if it is in the tmplist
                if (tmplist.includes(tmpnam)) {
                    console.log("General name of file " + files[i].name + " is in the processed list.");
                    continue;
                } else {
                    // search for the other occurrence
                    if (i == 0) {
                        console.log("File " + files[i].name + " is alone, though it seems to be a paired-end read.");
                        // There's no chance...
                        outlist.push([i]);
                        continue;
                    } else {
                        // There is chance!
                        let found = false;
                        for (let j = i - 1; j >= 0; j--) {
                            if (regExpWithTwoNumbers.test(files[j].name)) {
                                if (files[j].name === files[i].name) {
                                    console.log("File " + files[i].name + " has been uploaded at least twice. It will only be processed once.");
                                    tmplist.push(tmpnam); // To prevent multiple processings of same file.
                                } else {
                                    found = true;
                                    outlist.push([i, j]);
                                    tmplist.push(tmpnam);
                                    break;
                                }
                            }
                        }

                        if (!found) {
                            console.log("Compannion of file " + files[i].name + " not found. It will be processed alone.");
                            outlist.push([i]);
                        }
                    }
                }
            } else {
                console.log("File " + files[i].name + " is a single file read.");
                // Nope. It's alone!
                outlist.push([i]);
            }
        }

    }

    return outlist;
};
