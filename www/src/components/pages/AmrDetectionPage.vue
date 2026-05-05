<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <Pill class="w-6 h-6" />
        AMR detection
      </h1>

      <AmrDetectionHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1">
                <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Minimum fraction of stored gene-specific diagnostic k-mers required for an exact gene call.</p>
                </TooltipContent>
              </Tooltip>
              Minimum gene fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow" v-model="minGeneFractionPct" :lazy="true" :min="1" :max="25" :interval="1" :disabled="isDetectingAmr" />
              <span class="block w-[48px] text-center border border-gray-300 rounded-md text-sm">{{ minGeneFractionPct }}%</span>
            </div>
          </div>
          <div>
            <p class="flex items-center gap-1">
                <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Minimum fraction of stored family-specific diagnostic k-mers required for a family-level fallback call.</p>
                </TooltipContent>
              </Tooltip>
              Minimum family fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow" v-model="minFamilyFractionPct" :lazy="true" :min="5" :max="60" :interval="5" :disabled="isDetectingAmr" />
              <span class="block w-[48px] text-center border border-gray-300 rounded-md text-sm">{{ minFamilyFractionPct }}%</span>
            </div>
          </div>
        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 flex-1 overflow-hidden pt-12">
      <div v-if="amrError" class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="amrError === 'index'">Failed to load the AMR index file.</template>
        <template v-else-if="amrError === 'memory'">Error during processing, likely due to memory pressure.</template>
        <template v-else>An unexpected error occurred. Please reset and try again.</template>
      </div>

      <div v-if="!amrIndexLoaded && !isLoadingAmrIndex" v-bind="getRootPropsIndex()" :class="dropzoneClass">
        <input v-bind="getInputPropsIndex()" />
        <FileUp />
        <p v-if="isDragActiveIndex">Drop file here ...</p>
        <p v-else>Drop or click to upload your <b>AMR index file (.amridx)</b></p>
      </div>

      <div v-else-if="isLoadingAmrIndex" class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
        <Loader2 class="w-6 h-6 text-amber-500 animate-spin" />
        <p class="text-sm text-gray-500">Loading AMR index...</p>
      </div>

      <div v-else class="flex items-center gap-2 py-2 px-3 mx-6 bg-gray-50 rounded-md mb-4">
        <Check class="w-4 h-4 text-green-500" />
        <span class="flex-grow text-sm font-mono truncate">{{ allResults_amr.indexFileName }}</span>
        <span v-if="allResults_amr.indexInfo" class="text-xs text-gray-400 ml-2">{{ allResults_amr.indexInfo }}</span>
      </div>

      <template v-if="amrIndexLoaded">
        <div v-bind="getRootPropsFasta()" :class="dropzoneClass">
          <input v-bind="getInputPropsFasta()" />
          <FileUp />
          <p v-if="isDragActiveFasta">Drop FASTA files here ...</p>
          <p v-else>Drop or click to upload your <b>FASTA file(s)</b></p>
        </div>

        <div v-if="uploadedSampleNames.length > 0" class="mx-6 mt-4 max-h-48 overflow-y-auto">
          <div v-for="sampleName in uploadedSampleNames" :key="sampleName" class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="isFileDone(sampleName)" class="w-4 h-4 text-green-500" />
            <Loader2 v-else-if="isFileInFlight(sampleName)" class="w-4 h-4 text-amber-500 animate-spin" />
            <span class="flex-grow text-sm font-mono truncate">{{ sampleName }}</span>
          </div>
        </div>

        <div v-if="amrDetected" class="mx-6 mt-4 overflow-x-auto">
          <table class="w-full text-sm border border-gray-200 rounded-md">
            <thead>
              <tr class="bg-gray-50 text-left">
                <th class="px-3 py-2 font-medium text-gray-700">Sample</th>
                <th class="px-3 py-2 font-medium text-gray-700">Contig</th>
                <th class="px-3 py-2 font-medium text-gray-700">Call</th>
                <th class="px-3 py-2 font-medium text-gray-700">Class</th>
                <th class="px-3 py-2 font-medium text-gray-700">Region</th>
                <th class="px-3 py-2 font-medium text-gray-700">Hits</th>
                <th class="px-3 py-2 font-medium text-gray-700">Diag. Fraction</th>
                <th class="px-3 py-2 font-medium text-gray-700">Coverage</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="row in flattenedHits" :key="row.key" class="border-t border-gray-100">
                <td class="px-3 py-2 font-mono truncate max-w-[160px]">{{ row.sampleName }}</td>
                <td class="px-3 py-2 font-mono truncate max-w-[200px]">{{ row.hit.contig }}</td>
                <td class="px-3 py-2">{{ row.hit.gene_id ?? row.hit.gene_family }}</td>
                <td class="px-3 py-2">{{ row.hit.class_name ?? 'n/a' }}</td>
                <td class="px-3 py-2">{{ row.hit.start }}-{{ row.hit.end }}</td>
                <td class="px-3 py-2">{{ row.hit.distinct_hit_kmers }}/{{ row.hit.diagnostic_kmer_total }}</td>
                <td class="px-3 py-2">{{ (row.hit.diagnostic_kmer_fraction * 100).toFixed(1) }}%</td>
                <td class="px-3 py-2">{{ row.hit.reference_coverage_pct.toFixed(1) }}</td>
              </tr>
            </tbody>
          </table>
        </div>

        <div v-if="amrIndexLoaded" class="mx-6 mt-6">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
        </div>
      </template>
    </div>
  </div>
</template>

<script lang="ts">
import { computed, defineComponent, ref, Ref } from "vue";
import { useStore } from "vuex";
import { useDropzone } from "vue3-dropzone";
import { useActions, useGetters, useState } from "vuex-composition-helpers";
import { Check, FileUp, Info, Loader2, Pill, Trash2 } from "lucide-vue-next";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import VueSlider from "vue-3-slider-component";
import AmrDetectionHelpCollapsible from "@/components/help/AmrDetectionHelpCollapsible.vue";
import { Button } from "@/components/ui/button";
import { regExpForAnyFastx } from "@/utils";
import { AmrHit } from "@/types";

export default defineComponent({
  name: "AmrDetectionPage",
  components: {
    AmrDetectionHelpCollapsible,
    Button,
    Check,
    FileUp,
    Info,
    Loader2,
    Pill,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    Trash2,
    VueSlider,
  },
  props: {
    tabName: {
      type: String,
      required: true,
    },
  },
  setup() {
    useStore();
    const minGeneFractionPct: Ref<number> = ref(5);
    const minFamilyFractionPct: Ref<number> = ref(30);
    const uploadedSampleNames: Ref<string[]> = ref([]);

    const { loadAmrIndex, detectAmr, resetAllResults_amr } = useActions([
      "loadAmrIndex",
      "detectAmr",
      "resetAllResults_amr",
    ]);
    const { amrIndexLoaded, isLoadingAmrIndex, isDetectingAmr, detectingAmrFiles, amrDetected, amrResults, amrError } = useGetters([
      "amrIndexLoaded",
      "isLoadingAmrIndex",
      "isDetectingAmr",
      "detectingAmrFiles",
      "amrDetected",
      "amrResults",
      "amrError",
    ]) as any;
    const { allResults_amr } = useState(["allResults_amr"]) as any;

    const onDropIndex = (acceptFiles: File[]) => {
      if (acceptFiles.length > 0) {
        loadAmrIndex({ file: acceptFiles[0] });
      }
    };
    const onDropFasta = (acceptFiles: File[]) => {
      if (acceptFiles.length > 0) {
        uploadedSampleNames.value = acceptFiles.map((file) => file.name.replace(regExpForAnyFastx, ""));
        detectAmr({
          files: acceptFiles,
          min_gene_fraction: minGeneFractionPct.value / 100,
          min_family_fraction: minFamilyFractionPct.value / 100,
        });
      }
    };

    const { getRootProps: getRootPropsIndex, getInputProps: getInputPropsIndex, isDragActive: isDragActiveIndex } =
      useDropzone({ onDrop: onDropIndex, multiple: false });
    const { getRootProps: getRootPropsFasta, getInputProps: getInputPropsFasta, isDragActive: isDragActiveFasta } =
      useDropzone({ onDrop: onDropFasta, multiple: true });

    const dropzoneClass = [
      "p-6 mx-6 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600",
      "cursor-pointer hover:border-gray-400",
    ];

    const flattenedHits = computed(() => {
      const rows: { key: string; sampleName: string; hit: AmrHit }[] = [];
      const results = amrResults.value as Record<string, { hits: AmrHit[] }>;
      for (const [sampleName, result] of Object.entries(results)) {
        for (const hit of result.hits) {
          rows.push({ key: `${sampleName}-${hit.contig}-${hit.start}-${hit.gene_id ?? hit.gene_family}`, sampleName, hit });
        }
      }
      return rows;
    });

    const isFileDone = (sampleName: string) => sampleName in (allResults_amr.value.results || {});
    const isFileInFlight = (sampleName: string) => detectingAmrFiles.value.has(sampleName);
    const resetAll = () => {
      uploadedSampleNames.value = [];
      resetAllResults_amr();
    };

    return {
      allResults_amr,
      amrDetected,
      amrError,
      amrIndexLoaded,
      amrResults,
      detectingAmrFiles,
      dropzoneClass,
      flattenedHits,
      getInputPropsFasta,
      getInputPropsIndex,
      getRootPropsFasta,
      getRootPropsIndex,
      isDetectingAmr,
      isDragActiveFasta,
      isDragActiveIndex,
      isFileDone,
      isFileInFlight,
      isLoadingAmrIndex,
      minFamilyFractionPct,
      minGeneFractionPct,
      resetAll,
      uploadedSampleNames,
    };
  },
});
</script>
