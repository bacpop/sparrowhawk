<template>
  <div class="flex flex-row">
    <div class="w-1/3">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <Funnel class="w-6 h-6" />
        Host depletion
      </h1>

      <HostDepletionHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">

          <div class="flex flex-row items-center w-full gap-2">
            <input id="deplete" type="checkbox" v-model="deplete" :disabled="isFilteringDeacon"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">When enabled, reads matching the host index are removed. When disabled, matching reads are kept instead.</p>
              </TooltipContent>
            </Tooltip>
            <label for="deplete">
              Deplete mode
            </label>
          </div>

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Minimum number of matching k-mers required to classify a read as host-derived. Must be a positive integer.</p>
                </TooltipContent>
              </Tooltip>
              Absolute threshold
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="abs_threshold"
                         :lazy="true"
                         :min="1"
                         :max="50"
                         :interval="1"
                         :disabled="isFilteringDeacon"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ abs_threshold }}
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
                  <p class="max-w-xs">Minimum proportion of k-mers matching the host index required to classify a read as host-derived. A value between 0 and 1.</p>
                </TooltipContent>
              </Tooltip>
              Relative threshold
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="rel_threshold"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.01"
                         :disabled="isFilteringDeacon"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ rel_threshold }}
              </span>
            </div>
          </div>

        </div>
      </TooltipProvider>
    </div>

    <div class="w-2/3 pt-12">

      <div v-if="deaconError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="deaconError === 'memory'">
          Error during processing — most likely a memory issue. Try with fewer or smaller files.
        </template>
        <template v-else-if="deaconError === 'index'">
          Failed to load index file — the file may be corrupted or in the wrong format.
        </template>
        <template v-else>
          An unexpected error occurred. Please reset and try again.
        </template>
      </div>

      <div v-if="tabName === 'HostDepletion'">

        <!-- Phase 1: Index upload -->
        <div v-if="!deaconIndexLoaded && !isLoadingDeaconIndex"
             v-bind="getRootPropsIndex()"
             :class="[
               'p-6 mx-6 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind="getInputPropsIndex()" />
          <FileUp />
          <p v-if="isDragActiveIndex">Drop file here ...</p>
          <p v-else>Drop or click to upload your <b>Deacon index file (.idx)</b></p>
        </div>

        <div v-else-if="isLoadingDeaconIndex"
             class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
          <Loader2 class="w-6 h-6 text-amber-500 animate-spin" />
          <p class="text-sm text-gray-500">Loading index...</p>
        </div>

        <div v-else-if="deaconIndexLoaded"
             class="flex items-center gap-2 py-2 px-3 mx-6 bg-gray-50 rounded-md mb-4">
          <Check class="w-4 h-4 text-green-500" />
          <span class="flex-grow text-sm font-mono truncate">{{ allResults_deacon.indexFileName }}</span>
          <span v-if="allResults_deacon.indexInfo" class="text-xs text-gray-400 ml-2">{{ allResults_deacon.indexInfo }}</span>
        </div>

        <!-- Phase 2: Reads upload (only when index is loaded) -->
        <template v-if="deaconIndexLoaded && !isLoadingDeaconIndex">

          <div v-bind="getRootPropsReads()"
               :class="[
                 'p-6 mx-6 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
                 'cursor-pointer hover:border-gray-400'
               ]">
            <input v-bind="getInputPropsReads()" />
            <FileUp />
            <p v-if="isDragActiveReads">Drop files here ...</p>
            <p v-else>Drop or click to upload your <b>reads file(s) (.fastq)</b></p>
          </div>

          <!-- Per-sample status list -->
          <div v-if="uploadedSampleNames.length > 0" class="mx-6 mt-4 max-h-48 overflow-y-auto">
            <div v-for="sampleName in uploadedSampleNames" :key="sampleName"
                 class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
              <Check v-if="isFileDone(sampleName)" class="w-4 h-4 text-green-500"/>
              <Loader2 v-else-if="isFileInFlight(sampleName)" class="w-4 h-4 text-amber-500 animate-spin"/>
              <span class="flex-grow text-sm font-mono truncate">{{ sampleName }}</span>
            </div>
          </div>

          <!-- Bulk download (2+ results) -->
          <div v-if="resultCount >= 2" class="mx-6 mt-4 flex gap-2">
            <Button variant="outline" size="sm" @click="downloadZip">
              <Download class="mr-2 h-4 w-4" /> Download all (.zip)
            </Button>
            <Button variant="outline" size="sm" @click="downloadTarGz">
              <Download class="mr-2 h-4 w-4" /> Download all (.tar.gz)
            </Button>
          </div>

          <!-- Results table -->
          <div v-if="deaconFiltered" class="mx-6 mt-4 overflow-x-auto">
            <table class="w-full text-sm border border-gray-200 rounded-md">
              <thead>
                <tr class="bg-gray-50 text-left">
                  <th class="px-3 py-2 font-medium text-gray-700">Sample name</th>
                  <th class="px-3 py-2 font-medium text-gray-700">Total reads</th>
                  <th class="px-3 py-2 font-medium text-gray-700">Kept reads</th>
                  <th class="px-3 py-2 font-medium text-gray-700">Removed reads</th>
                  <th class="px-3 py-2 font-medium text-gray-700">Download</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="(result, sampleName) in deaconResults" :key="sampleName"
                    class="border-t border-gray-100">
                  <td class="px-3 py-2 font-mono truncate max-w-[200px]">{{ result.sampleName }}</td>
                  <td class="px-3 py-2">{{ result.totalReads }}</td>
                  <td class="px-3 py-2">{{ result.keptReads }}</td>
                  <td class="px-3 py-2">{{ result.removedReads }}</td>
                  <td class="px-3 py-2 flex gap-1 flex-wrap">
                    <Button variant="outline" size="sm" @click="downloadR1(result)">
                      <Download class="mr-1 h-3 w-3" />
                      {{ result.outputGzip2 ? 'R1 ' : '' }}(.fastq.gz)
                    </Button>
                    <Button v-if="result.outputGzip2" variant="outline" size="sm" @click="downloadR2(result)">
                      <Download class="mr-1 h-3 w-3" /> R2 (.fastq.gz)
                    </Button>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>

        </template>

        <!-- Reset button -->
        <div v-if="deaconIndexLoaded" class="mx-6 mt-6">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
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
import { Check, FileUp, Loader2, Funnel, Info, Download, Trash2 } from "lucide-vue-next";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import VueSlider from 'vue-3-slider-component';
import HostDepletionHelpCollapsible from "@/components/help/HostDepletionHelpCollapsible.vue";
import { Button } from "@/components/ui/button";
import { fastqExtensionsWithDotAndCompressList, getFilesToProcess, regExpWithTwoNumbers, regExpForAnyFastx } from "@/utils";
import { buildZipBinary, buildTarGzBinary } from "@/archiveUtils";
import { DepletionResult } from "@/types";

export default defineComponent({
  name: "HostDepletionPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    Check,
    FileUp,
    Loader2,
    Funnel,
    Info,
    Download,
    Trash2,
    VueSlider,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    HostDepletionHelpCollapsible,
    Button,
  },
  setup() {
    const store = useStore();

    const deplete: Ref<boolean> = ref(true);
    const abs_threshold: Ref<number> = ref(1);
    const rel_threshold: Ref<number> = ref(0.05);
    const uploadedSampleNames: Ref<string[]> = ref([]);

    const { loadDeaconIndex, filterDeaconReads, resetAllResults_deacon } = useActions([
      "loadDeaconIndex",
      "filterDeaconReads",
      "resetAllResults_deacon",
    ]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const { allResults_deacon } = useState(["allResults_deacon"]) as any;

    const deaconResults = computed<Record<string, DepletionResult>>(() => store.getters.deaconResults);
    const resultCount = computed<number>(() => Object.keys(deaconResults.value).length);

    function resetAll(): void {
      uploadedSampleNames.value = [];
      resetAllResults_deacon();
    }

    function onDropIndex(acceptFiles: File[]): void {
      if (acceptFiles.length === 0) return;
      loadDeaconIndex({ file: acceptFiles[0] });
    }

    function onDropReads(acceptFiles: File[]): void {
      if (acceptFiles.length === 0) return;
      // Accumulate sample names from this batch
      const indxlist = getFilesToProcess(acceptFiles);
      for (const sublist of indxlist) {
        const file = acceptFiles[sublist[0]];
        const sampleName = sublist.length > 1
          ? file.name.replace(regExpWithTwoNumbers, '')
          : file.name.replace(regExpForAnyFastx, '');
        if (!uploadedSampleNames.value.includes(sampleName)) {
          uploadedSampleNames.value.push(sampleName);
        }
      }
      filterDeaconReads({
        files: acceptFiles,
        deplete: deplete.value,
        abs_threshold: abs_threshold.value,
        rel_threshold: rel_threshold.value,
      });
    }

    function downloadR1(result: DepletionResult): void {
      const name = result.sampleName + (result.outputGzip2 ? '_R1' : '') + '_filtered.fastq.gz';
      const blob = new Blob([result.outputGzip.buffer as ArrayBuffer], { type: 'application/gzip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = name; a.click();
      URL.revokeObjectURL(url);
    }

    function downloadR2(result: DepletionResult): void {
      if (!result.outputGzip2) return;
      const name = result.sampleName + '_R2_filtered.fastq.gz';
      const blob = new Blob([result.outputGzip2.buffer as ArrayBuffer], { type: 'application/gzip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = name; a.click();
      URL.revokeObjectURL(url);
    }

    function downloadZip(): void {
      const files: Record<string, Uint8Array> = {};
      for (const result of Object.values(deaconResults.value)) {
        files[result.sampleName + (result.outputGzip2 ? '_R1' : '') + '_filtered.fastq.gz'] = result.outputGzip;
        if (result.outputGzip2) {
          files[result.sampleName + '_R2_filtered.fastq.gz'] = result.outputGzip2;
        }
      }
      const bytes = buildZipBinary(files);
      const blob = new Blob([bytes.buffer as ArrayBuffer], { type: 'application/zip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = 'depleted.zip'; a.click();
      URL.revokeObjectURL(url);
    }

    function downloadTarGz(): void {
      const files: Record<string, Uint8Array> = {};
      for (const result of Object.values(deaconResults.value)) {
        files[result.sampleName + (result.outputGzip2 ? '_R1' : '') + '_filtered.fastq.gz'] = result.outputGzip;
        if (result.outputGzip2) {
          files[result.sampleName + '_R2_filtered.fastq.gz'] = result.outputGzip2;
        }
      }
      const bytes = buildTarGzBinary(files);
      const blob = new Blob([bytes.buffer as ArrayBuffer], { type: 'application/gzip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = 'depleted.tar.gz'; a.click();
      URL.revokeObjectURL(url);
    }

    const {
      getRootProps: getRootPropsIndex,
      getInputProps: getInputPropsIndex,
      isDragActive: isDragActiveIndex,
    } = useDropzone({
      onDrop: onDropIndex,
      accept: [".idx"],
      multiple: false,
    });

    const {
      getRootProps: getRootPropsReads,
      getInputProps: getInputPropsReads,
      isDragActive: isDragActiveReads,
    } = useDropzone({
      onDrop: onDropReads,
      accept: fastqExtensionsWithDotAndCompressList,
      multiple: true,
    });

    return {
      deplete,
      abs_threshold,
      rel_threshold,
      allResults_deacon,
      uploadedSampleNames,
      deaconResults,
      resultCount,
      store,
      resetAll,
      downloadR1,
      downloadR2,
      downloadZip,
      downloadTarGz,
      getRootPropsIndex,
      getInputPropsIndex,
      isDragActiveIndex,
      getRootPropsReads,
      getInputPropsReads,
      isDragActiveReads,
    };
  },
  computed: {
    deaconError(): string | null {
      return this.store.getters.deaconError;
    },
    deaconIndexLoaded(): boolean {
      return this.store.getters.deaconIndexLoaded;
    },
    isLoadingDeaconIndex(): boolean {
      return this.store.getters.isLoadingDeaconIndex;
    },
    isFilteringDeacon(): boolean {
      return this.store.getters.isFilteringDeacon;
    },
    deaconFiltered(): boolean {
      return this.store.getters.deaconFiltered;
    },
    filteringDeaconFilesSet(): Set<string> {
      return this.store.getters.filteringDeaconFiles;
    },
  },
  methods: {
    isFileDone(sampleName: string): boolean {
      return sampleName in (this.deaconResults || {});
    },
    isFileInFlight(sampleName: string): boolean {
      return this.filteringDeaconFilesSet.has(sampleName);
    },
  },
});
</script>

<style scoped>
</style>
