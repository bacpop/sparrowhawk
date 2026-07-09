<template>
  <div class="flex flex-col gap-6 md:flex-row md:gap-0">
    <div class="w-full md:w-[350px] md:shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <ChartNetwork class="w-6 h-6" />
        Transmission
      </h1>

      <TransmissionHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1 mb-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Maximum pairwise SNP distance for two samples to be placed in the same transmission cluster.</p>
                </TooltipContent>
              </Tooltip>
              SNP threshold
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider
                class="flex-grow"
                v-model="snp_threshold"
                :lazy="true"
                :min="1"
                :max="1000"
                :interval="1"
                :disabled="isTransmissionStandaloneClustering"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ snp_threshold }}
              </span>
            </div>
          </div>

          <Button
            :disabled="!alignmentFile || isTransmissionStandaloneClustering"
            @click="runClustering"
            class="w-full"
            variant="outline"
            size="sm"
          >
            <Network class="mr-1 h-3 w-3" />
            Run transmission clustering
          </Button>

          <hr class="border-gray-200" />

          <div>
            <p class="text-sm mb-1">Optional metadata file (CSV)</p>
            <label class="block cursor-pointer">
              <span class="block p-2 border border-gray-200 rounded-md text-sm text-gray-500 hover:border-gray-400 text-center">
                {{ metadataFileName || 'Upload metadata CSV' }}
              </span>
              <input type="file" accept=".csv" class="hidden" @change="onMetadataFileChange" />
            </label>
            <div v-if="metadataError" class="mt-2 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
              {{ metadataError }}
            </div>
          </div>

          <div v-if="metadataRows.length > 0" class="flex items-center gap-2">
            <input id="transmission_enable_location" type="checkbox" v-model="enableLocationMatching" />
            <label for="transmission_enable_location" class="text-sm">Enable location matching</label>
          </div>
        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 w-full flex-1 overflow-visible pt-0 md:overflow-hidden md:pt-12">
      <div v-if="transmissionStandaloneError" class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        {{ transmissionErrorMessage }}
      </div>

      <div
        v-if="!isTransmissionStandaloneClustering"
        v-bind="getRootPropsTransmission()"
        :class="[
          'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
          'cursor-pointer hover:border-gray-400'
        ]"
      >
        <input v-bind="getInputPropsTransmission()" />
        <FileUp />
        <p v-if="isDragActiveTransmission">Drop alignment here ...</p>
        <p v-else>Drop or click to upload an <b>aligned FASTA/ALN</b></p>
      </div>

      <div v-else class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
        <Loader2 class="w-6 h-6 text-amber-500 animate-spin" />
        <p class="text-sm text-gray-500">Clustering...</p>
      </div>

      <div v-if="alignmentFileName" class="mx-6 mr-0 mt-4 max-h-48 overflow-y-auto">
        <div class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
          <Loader2 v-if="isTransmissionStandaloneClustering" class="w-4 h-4 text-orange-500 animate-spin" />
          <Check v-else-if="hasTransmissionStandaloneClusterResults" class="w-4 h-4 text-green-500" />
          <span class="flex-grow text-sm font-mono truncate">{{ alignmentFileName }}</span>
        </div>
      </div>

      <div v-if="alignmentFile" class="mx-6 mr-0 mt-4 flex items-center gap-2">
        <Button @click="clear" variant="outline" size="sm">
          <Trash2 class="mr-1 h-3 w-3" />
          Clear results
        </Button>
      </div>

      <TransmissionClusterResults
        v-if="hasTransmissionStandaloneClusterResults"
        class="mx-6 mt-6"
        :clusterResults="transmissionStandalone.clusterResults"
        :transmissionGraph="transmissionStandalone.transmissionGraph"
        :metadataRows="metadataRows"
        :enableLocationMatching="enableLocationMatching"
      />
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, Ref } from "vue";
import { useDropzone } from "vue3-dropzone";
import { useActions, useState } from "vuex-composition-helpers";
import { useStore } from "vuex";
import VueSlider from "vue-3-slider-component";
import { ChartNetwork, Check, FileUp, Info, Loader2, Network, Trash2 } from "lucide-vue-next";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import TransmissionClusterResults from "@/components/TransmissionClusterResults.vue";
import TransmissionHelpCollapsible from "@/components/help/TransmissionHelpCollapsible.vue";
import { parseTransmissionMetadataCsv } from "@/utils/transmissionMetadata";
import type { MetadataRow } from "@/types";

const alignmentExtensions = [".fa", ".fasta", ".fas", ".fna", ".aln"];

export default defineComponent({
  name: "TransmissionPage",
  props: {
    tabName: {
      type: String,
      required: true,
    },
  },
  components: {
    Button,
    ChartNetwork,
    Check,
    FileUp,
    Info,
    Loader2,
    Network,
    Trash2,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    TransmissionClusterResults,
    TransmissionHelpCollapsible,
    VueSlider,
  },
  setup() {
    const store = useStore();
    const snp_threshold: Ref<number> = ref(20);
    const alignmentFile: Ref<File | null> = ref(null);
    const alignmentFileName: Ref<string | null> = ref(null);
    const metadataRows: Ref<MetadataRow[]> = ref([]);
    const metadataError: Ref<string | null> = ref(null);
    const metadataFileName: Ref<string | null> = ref(null);
    const enableLocationMatching: Ref<boolean> = ref(false);

    const { processTransmissionStandaloneCluster, resetTransmissionStandaloneResults } = useActions([
      "processTransmissionStandaloneCluster",
      "resetTransmissionStandaloneResults",
    ]);
    const { transmissionStandalone } = useState(["transmissionStandalone"]) as any;

    function onDropTransmission(acceptFiles: File[]): void {
      const file = acceptFiles[0];
      if (!file) {
        return;
      }
      alignmentFile.value = file;
      alignmentFileName.value = file.name;
      resetTransmissionStandaloneResults();
    }

    function parseMetadataCSV(file: File): void {
      metadataError.value = null;
      metadataFileName.value = file.name;
      const reader = new FileReader();
      reader.onload = (event) => {
        try {
          metadataRows.value = parseTransmissionMetadataCsv(String(event.target?.result ?? ""));
        } catch (error) {
          metadataRows.value = [];
          metadataError.value = error instanceof Error ? error.message : String(error);
        }
      };
      reader.readAsText(file);
    }

    function onMetadataFileChange(event: Event): void {
      const input = event.target as HTMLInputElement;
      if (input.files && input.files[0]) {
        parseMetadataCSV(input.files[0]);
      }
    }

    function clear(): void {
      alignmentFile.value = null;
      alignmentFileName.value = null;
      snp_threshold.value = 20;
      metadataRows.value = [];
      metadataError.value = null;
      metadataFileName.value = null;
      enableLocationMatching.value = false;
      resetTransmissionStandaloneResults();
    }

    const {
      getRootProps: getRootPropsTransmission,
      getInputProps: getInputPropsTransmission,
      isDragActive: isDragActiveTransmission,
      ...restTransmission
    } = useDropzone({
      onDrop: onDropTransmission,
      accept: alignmentExtensions,
      multiple: false,
    });

    return {
      store,
      snp_threshold,
      alignmentFile,
      alignmentFileName,
      metadataRows,
      metadataError,
      metadataFileName,
      enableLocationMatching,
      processTransmissionStandaloneCluster,
      transmissionStandalone,
      getRootPropsTransmission,
      getInputPropsTransmission,
      isDragActiveTransmission,
      onMetadataFileChange,
      clear,
      ...restTransmission,
    };
  },
  computed: {
    isTransmissionStandaloneClustering(): boolean {
      return this.store.getters.isTransmissionStandaloneClustering;
    },
    hasTransmissionStandaloneClusterResults(): boolean {
      return this.store.getters.hasTransmissionStandaloneClusterResults;
    },
    transmissionStandaloneError(): string | null {
      return this.store.getters.transmissionStandaloneError;
    },
    transmissionErrorMessage(): string {
      if (this.transmissionStandaloneError === 'worker_unavailable') {
        return 'Transmission clustering worker is not available. Please reload and try again.';
      }
      if (this.transmissionStandaloneError === 'memory') {
        return 'Error during processing, most likely a memory issue. Try with a smaller alignment.';
      }
      return this.transmissionStandaloneError || 'An unexpected error occurred. Please reset and try again.';
    },
  },
  methods: {
    runClustering(): void {
      if (!this.alignmentFile) {
        return;
      }
      this.processTransmissionStandaloneCluster({
        file: this.alignmentFile,
        snp_threshold: this.snp_threshold,
      });
    },
  },
});
</script>
