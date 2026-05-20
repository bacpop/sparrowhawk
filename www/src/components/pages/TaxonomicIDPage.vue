<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <ScanFace class="w-6 h-6" />
        Taxonomic ID
      </h1>

      <TaxonomicIDHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Controls the filtering of nucleotides depending on the sequencing error information.</p>
                </TooltipContent>
              </Tooltip>
              Min Illumina read quality
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_qual"
                         :lazy="true"
                         :min="0"
                         :max="33"
                         :interval="1"
                         :disabled="isIdentifying"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_qual }}
              </span>
            </div>
          </div>

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Only k-mers appearing more than this value will be used.</p>
                </TooltipContent>
              </Tooltip>
              Min counts for k-mer filtering
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_count"
                         :lazy="true"
                         :min="1"
                         :max="30"
                         :interval="1"
                         :disabled="isIdentifying"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_count }}
              </span>
            </div>
          </div>

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">A real number between 0 and 1 that controls what proportion of reads to use for processing.</p>
                </TooltipContent>
              </Tooltip>
              Proportion of reads
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="proportion_reads"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.05"
                         :disabled="isIdentifying"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ proportion_reads }}
              </span>
            </div>
          </div>

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Number of parallel workers used to process samples simultaneously. Higher values speed up processing of multiple files but use more memory.</p>
                </TooltipContent>
              </Tooltip>
              Workers
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="numWorkers"
                         :lazy="true"
                         :min="1"
                         :max="8"
                         :interval="1"
                         :disabled="isIdentifying"
                         @change="onNumWorkersChange"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ numWorkers }}
              </span>
            </div>
          </div>
        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 flex-1 overflow-hidden pt-12">

      <div v-if="sketchlibError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="sketchlibError === 'memory'">
          Error during processing — most likely a memory issue. Try with fewer or smaller files.
        </template>
        <template v-else>
          An unexpected error occurred. Please reset and try again.
        </template>
      </div>

      <div v-if="tabName=='TaxonomicID'">

        <!-- Upload dropbox - always visible when not identifying -->
        <div v-if="!isIdentifying"
             v-bind='getRootPropsSample()'
             :class="[
               'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind='getInputPropsSample()' />
          <FileUp/>
          <p v-if='isDragActiveSample'>
            Drop files here ...
          </p>
          <p v-else>
            Drop or click to upload your <b>sample fastq/a file(s)</b>
          </p>
        </div>

        <div v-else class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
          <Loader2 class="w-6 h-6 text-amber-500 animate-spin"/>
          <p class="text-sm text-gray-500">
            Identifying {{ identifyingFilesArray.length }} sample(s)...
          </p>
          <div v-if="identifyingFilesArray.length > 0" class="text-xs text-gray-400 mt-1">
            Processing: {{ identifyingFilesArray.join(', ') }}
          </div>
        </div>

        <!-- Show uploaded files with per-file status -->
        <div v-if="uploadedFileNames.length > 0" class="mx-6 mr-0 mt-4 max-h-48 overflow-y-auto">
          <div v-for="fileName in uploadedFileNames" :key="fileName"
               class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="isFileIdentified(fileName)" class="w-4 h-4 text-green-500"/>
            <Loader2 v-else-if="isFileIdentifying(fileName)" class="w-4 h-4 text-amber-500 animate-spin"/>
            <span class="flex-grow text-sm font-mono truncate">
              {{ fileName }}
            </span>
          </div>
        </div>

        <!-- Clear button -->
        <div v-if="uploadedFileNames.length > 0" class="mx-6 mr-0 mt-4">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
        </div>

        <!-- Show results as a single table with expandable rows per sample -->
        <div v-if="sampleIdentified" class="px-6 mt-4">
          <div class="flex justify-start mb-2">
            <Button variant="outline" size="sm" @click="downloadTsv">
              <Download class="mr-2 h-4 w-4" />
              Download TSV
            </Button>
          </div>
          <div class="max-h-96 overflow-y-auto">
            <DataTable :columns="tableColumns" :data="tableData" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, Ref, computed } from "vue";
import { useStore } from "vuex";
import { useDropzone } from "vue3-dropzone";
import { useActions, useState } from "vuex-composition-helpers";
import VueSlider from 'vue-3-slider-component';
import { Check, FileUp, Loader2, Info, ScanFace, Download, Trash2 } from "lucide-vue-next";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import TaxonomicIDHelpCollapsible from "@/components/help/TaxonomicIDHelpCollapsible.vue";
import DataTable from "@/components/pages/taxonomic-id/DataTable.vue";
import { columns, TaxonomicIDRow } from "@/components/pages/taxonomic-id/columns";
import { SampleIdentifyResult } from "@/types";
import { fastxExtensionsWithDotAndCompressList } from "@/utils";

export default defineComponent({
  name: "TaxonomicIDPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    VueSlider,
    Check,
    FileUp,
    Loader2,
    Info,
    ScanFace,
    Download,
    Trash2,
    Button,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    TaxonomicIDHelpCollapsible,
    DataTable
  },
  setup() {
    const store = useStore();
    const proportion_reads: Ref<number> = ref(1);
    const min_qual: Ref<number> = ref(20);
    const min_count: Ref<number> = ref(5);
    const uploadedFileNames: Ref<string[]> = ref([]);

    const { identifyFiles, resetAllResults_sketchlib, initSketchlibWorkers } = useActions(["identifyFiles", "resetAllResults_sketchlib", "initSketchlibWorkers"]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const { allResults_sketchlib } = useState(["allResults_sketchlib"]) as any;

    const numWorkers: Ref<number> = ref(4);

    function onNumWorkersChange(value: number): void {
      initSketchlibWorkers(value);
    }

    function onDropSample(acceptFiles: File[]): void {
      uploadedFileNames.value = acceptFiles.map(f => f.name);
      identifyFiles({ acceptFiles: acceptFiles, proportion_reads: proportion_reads.value, min_count: min_count.value, min_qual: min_qual.value });
    }
    function resetAll(): void {
      uploadedFileNames.value = [];
      resetAllResults_sketchlib();
    }

    // Get sample name from file name (same logic as action)
    //todo: need more testing
    function getSampleNameFromFile(fileName: string): string {
      if (/(_1|_2)(.fastq.gz|.fq.gz)$/.test(fileName)) {
        return fileName.replace(/(_1.fastq.gz|_1.fq.gz|_2.fastq.gz|_2.fq.gz)$/, '');
      }
      return fileName.replace(/(.fasta|.fasta.gz|.fna|.fna.gz|.fa|.fa.gz|.fq|.fq.gz|.fastq|.fastq.gz)$/, '');
    }

    const tableColumns = columns;

    const tableData = computed<TaxonomicIDRow[]>(() => {
      const results = allResults_sketchlib.value.results as Record<string, SampleIdentifyResult>;
      const topLevelRows: TaxonomicIDRow[] = [];
      for (const [sampleName, sampleResult] of Object.entries(results || {})) {
        if (!sampleResult?.idSpecies?.length) continue;

        const allRows: TaxonomicIDRow[] = sampleResult.idSpecies.map((species, i) => {
          const parts = (sampleResult.idMetadata[i] ?? '').split('|');
          return {
            sample: sampleName,
            rank: i + 1,
            species,
            probability: sampleResult.idProbs[i],
            metaSpecies: parts[0]?.trim() ?? '',
            metaGemsparcl: parts[1]?.trim() ?? '',
            metaGtdb: parts[2]?.trim() ?? '',
          };
        });

        topLevelRows.push({
          ...allRows[0],
          clusterDetails: allRows,
        });
      }
      return topLevelRows;
    });

    function downloadTsv(): void {
      const headers = [
        "Sample", "Rank", "Species", "Similarity (%)",
        "Species (metadata)", "Gemsparcl ID", "GTDB species composition",
      ];

      const rows: string[][] = [];
      for (const row of tableData.value) {
        for (const sub of row.clusterDetails ?? [row]) {
          rows.push([
            sub.sample,
            String(sub.rank),
            sub.species,
            (sub.probability * 100).toFixed(2),
            sub.metaSpecies,
            sub.metaGemsparcl,
            sub.metaGtdb,
          ]);
        }
      }

      const tsv = [headers, ...rows].map(r => r.join("\t")).join("\n");
      const blob = new Blob([tsv], { type: "text/tab-separated-values" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "sparrowhawk_id_results.tsv";
      a.click();
      URL.revokeObjectURL(url);
    }

    const {
      getRootProps: getRootPropsSample,
      getInputProps: getInputPropsSample,
      isDragActive: isDragActiveSample,
      ...restSample
    } = useDropzone({
      onDrop: onDropSample,
      accept: fastxExtensionsWithDotAndCompressList,
      multiple: true
    });

    return {
      proportion_reads,
      min_qual,
      min_count,
      numWorkers,
      onNumWorkersChange,
      uploadedFileNames,
      resetAll,
      getRootPropsSample,
      getInputPropsSample,
      isDragActiveSample,
      onDropSample,
      getSampleNameFromFile,
      tableColumns,
      tableData,
      downloadTsv,
      allResults_sketchlib,
      store,
      ...restSample,
    };
  },
  computed: {
    sketchlibError(): string | null {
      return this.store.getters.sketchlibError;
    },
    sampleIdentified(): boolean {
      return this.store.getters.sampleIdentified;
    },
    isIdentifying(): boolean {
      return this.store.getters.isIdentifying;
    },
    identifyingFilesSet(): Set<string> {
      return this.store.getters.isIdentifyingFiles;
    },
    identifyingFilesArray(): string[] {
      return Array.from(this.identifyingFilesSet);
    }
  },

  methods: {
    clear(): void {
      this.resetAll();
    },
    isFileIdentified(fileName: string): boolean {
      const sampleName = this.getSampleNameFromFile(fileName);
      return sampleName in (this.allResults_sketchlib.results || {});
    },
    isFileIdentifying(fileName: string): boolean {
      const sampleName = this.getSampleNameFromFile(fileName);
      return this.identifyingFilesSet.has(sampleName);
    }
  },
});
</script>

<style scoped>
</style>
