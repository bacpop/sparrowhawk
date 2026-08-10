<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <Pill class="w-6 h-6" />
        AMR detection
      </h1>

      <AMRDetectionHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Minimum fraction of exact-gene diagnostic k-mers required to call a gene report unit.</p>
                </TooltipContent>
              </Tooltip>
              Minimum exact-gene fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_gene_fraction"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.01"
                         :disabled="isDetectingAmr"
              />
              <span class="block w-[46px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_gene_fraction.toFixed(2) }}
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
                  <p class="max-w-xs">Minimum fraction of diagnostic k-mers required to call collapsed gene-group or hierarchy report units.</p>
                </TooltipContent>
              </Tooltip>
              Minimum gene-group fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_gene_group_fraction"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.01"
                         :disabled="isDetectingAmr"
              />
              <span class="block w-[46px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_gene_group_fraction.toFixed(2) }}
              </span>
            </div>
          </div>
        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 flex-1 overflow-hidden pt-12">
      <div v-if="amrError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="amrError === 'file_count'">
          Upload exactly one FASTA file.
        </template>
        <template v-else-if="amrError === 'fasta'">
          AMR detection currently accepts FASTA assemblies only.
        </template>
        <template v-else-if="amrError === 'index'">
          Could not load the bundled AMR index.
        </template>
        <template v-else-if="amrError === 'worker'">
          AMR worker is not available yet. Please try again after the page finishes loading.
        </template>
        <template v-else>
          AMR detection failed: {{ amrError }}
        </template>
      </div>

      <div v-if="tabName=='AMRDetection'">
        <div v-if="!isDetectingAmr"
             v-bind='getRootPropsSample()'
             :class="[
               'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind='getInputPropsSample()' />
          <FileUp/>
          <p v-if='isDragActiveSample'>
            Drop file here ...
          </p>
          <p v-else>
            Drop or click to upload one <b>assembly FASTA file</b>
          </p>
        </div>

        <div v-else class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
          <Loader2 class="w-6 h-6 text-amber-500 animate-spin"/>
          <p class="text-sm text-gray-500">
            Detecting AMR...
          </p>
          <div v-if="detectingAmrFilesArray.length > 0" class="text-xs text-gray-400 mt-1">
            Processing: {{ detectingAmrFilesArray.join(', ') }}
          </div>
        </div>

        <div v-if="uploadedFileNames.length > 0" class="mx-6 mr-0 mt-4 max-h-48 overflow-y-auto">
          <div v-for="fileName in uploadedFileNames" :key="fileName"
               class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="amrDetected && !isFileDetecting(fileName)" class="w-4 h-4 text-green-500"/>
            <Loader2 v-else-if="isFileDetecting(fileName)" class="w-4 h-4 text-amber-500 animate-spin"/>
            <span class="flex-grow text-sm font-mono truncate">
              {{ fileName }}
            </span>
          </div>
        </div>

        <div v-if="uploadedFileNames.length > 0 || amrDetected" class="mx-6 mr-0 mt-4">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
        </div>

        <div v-if="amrDetected" class="px-6 mt-4">
          <div class="flex items-center justify-between mb-2 gap-3">
            <div class="text-sm text-gray-500 truncate">
              <span v-if="allResults_amr.indexFileName" class="font-mono">{{ allResults_amr.indexFileName }}</span>
              <span v-if="amrResult"> · {{ amrResult.hits.length }} hit(s)</span>
              <span v-if="amrElapsedLabel"> · {{ amrElapsedLabel }}</span>
              <span v-if="amrMemoryLabel"> · peak WebAssembly memory {{ amrMemoryLabel }}</span>
            </div>
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
import { Check, Download, FileUp, Info, Loader2, Pill, Trash2 } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import AMRDetectionHelpCollapsible from "@/components/help/AMRDetectionHelpCollapsible.vue";
import DataTable from "@/components/pages/taxonomic-id/DataTable.vue";
import { columns, AmrDetectionRow } from "@/components/pages/amr-detection/columns";
import { AmrDetectionResult } from "@/types";
import { buildAmrTsv } from "@/amrTsv";
import { fastaExtensionsWithDotAndCompressList, formatBytes, formatDuration, regExpForAnyFasta } from "@/utils";
import {saveTextFile} from "@/platform/files";

export default defineComponent({
  name: "AMRDetectionPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    VueSlider,
    Check,
    Download,
    FileUp,
    Info,
    Loader2,
    Pill,
    Trash2,
    Button,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    AMRDetectionHelpCollapsible,
    DataTable,
  },
  setup() {
    const store = useStore();
    const min_gene_fraction: Ref<number> = ref(0.10);
    const min_gene_group_fraction: Ref<number> = ref(0.10);
    const uploadedFileNames: Ref<string[]> = ref([]);

    const { detectAmrFile, resetAllResults_amr } = useActions(["detectAmrFile", "resetAllResults_amr"]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const { allResults_amr } = useState(["allResults_amr"]) as any;

    function onDropSample(acceptFiles: File[]): void {
      uploadedFileNames.value = acceptFiles.map(f => f.name);
      detectAmrFile({
        files: acceptFiles,
        min_gene_fraction: min_gene_fraction.value,
        min_gene_group_fraction: min_gene_group_fraction.value,
      });
    }

    function resetAll(): void {
      uploadedFileNames.value = [];
      resetAllResults_amr();
    }

    function getSampleNameFromFile(fileName: string): string {
      return fileName.replace(regExpForAnyFasta, '');
    }

    const tableColumns = columns;

    const amrResult = computed<AmrDetectionResult | null>(() => allResults_amr.value.result as AmrDetectionResult | null);

    const tableData = computed<AmrDetectionRow[]>(() => {
      const result = amrResult.value;
      if (!result) return [];
      return result.hits.map((hit) => ({
        sample: result.sample_name,
        query: hit.query_id,
        unitLabel: hit.unit_label,
        callType: hit.call_type,
        category: hit.type_name ?? '',
        subtype: hit.subtype ?? '',
        className: hit.class_name ?? '',
        subclass: hit.subclass ?? '',
        callFraction: hit.call_fraction,
        diagnosticKmers: `${hit.first_pass_distinct}/${hit.first_pass_diagnostic_total}`,
        span: `${hit.start}-${hit.end}`,
        memberCount: hit.member_count,
      }));
    });

    async function downloadTsv(): Promise<void> {
      const result = amrResult.value;
      if (!result) return;
      const tsv = buildAmrTsv(result);
      await saveTextFile(tsv, "sparrowhawk_amr_results.tsv", "text/tab-separated-values;charset=utf-8");
    }

    const {
      getRootProps: getRootPropsSample,
      getInputProps: getInputPropsSample,
      isDragActive: isDragActiveSample,
      ...restSample
    } = useDropzone({
      onDrop: onDropSample,
      accept: fastaExtensionsWithDotAndCompressList,
      multiple: false
    });

    return {
      min_gene_fraction,
      min_gene_group_fraction,
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
      allResults_amr,
      amrResult,
      store,
      ...restSample,
    };
  },
  computed: {
    amrError(): string | null {
      return this.store.getters.amrError;
    },
    amrDetected(): boolean {
      return this.store.getters.amrDetected;
    },
    isDetectingAmr(): boolean {
      return this.store.getters.isDetectingAmr;
    },
    detectingAmrFilesSet(): Set<string> {
      return this.store.getters.detectingAmrFiles;
    },
    detectingAmrFilesArray(): string[] {
      return Array.from(this.detectingAmrFilesSet);
    },
    amrElapsedLabel(): string {
      const ms = this.amrResult?.elapsedMs;
      return ms == null ? "" : formatDuration(ms);
    },
    amrMemoryLabel(): string {
      const bytes = this.amrResult?.wasmMemoryBytes;
      return bytes == null ? "" : formatBytes(bytes);
    }
  },
  methods: {
    isFileDetecting(fileName: string): boolean {
      const sampleName = this.getSampleNameFromFile(fileName);
      return this.detectingAmrFilesSet.has(sampleName);
    }
  },
});
</script>
