<template>
  <div class="flex flex-col gap-6 md:flex-row md:gap-0">
    <div class="w-full md:w-[350px] md:shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <ScanBox class="w-6 h-6" />
        Embed proteins
      </h1>

      <ProteinEmbeddingsHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div class="flex flex-row items-center w-full gap-2">
            <input id="use_gpu" type="checkbox" v-model="use_gpu" :disabled="isEmbedding || !gpuSupported"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Run the model on your GPU via WebGPU. This might accelerate the results, depending on browser, operating system, and hardware. Off by default; results are identical either way.</p>
              </TooltipContent>
            </Tooltip>
            <label for="use_gpu">
              GPU acceleration
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2" :class="!use_gpu ? 'opacity-50' : ''">
            <input id="auto_device" type="checkbox" v-model="auto_device" :disabled="!use_gpu || isEmbedding"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Let the browser pick. Uncheck to choose between the GPUs it reports. Useful when a dedicated/discrete card and an integrated one are both present.</p>
              </TooltipContent>
            </Tooltip>
            <label for="auto_device">
              Automatic device selection
            </label>
          </div>

          <div v-if="use_gpu && !auto_device">
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">The page cannot enumerate GPUs, but it can ask for a high-performance or a low-power one, and the browser honours that hint.</p>
                </TooltipContent>
              </Tooltip>
              GPU device
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSelect class="flex-grow"
                         v-model="gpu_power_pref"
                         :options="adapterOptions"
                         :is-disabled="isEmbedding"
              />
            </div>
            <p v-if="rendererHint" class="text-xs text-gray-400 mt-1">
              Your browser reports {{ rendererHint }} for graphics, but will not say which of
              these devices that is.
            </p>
          </div>

          <div v-if="use_gpu">
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">How many pieces of GPU work are queued into a single
                    submission. Changing this reloads the model.</p>
                </TooltipContent>
              </Tooltip>
              GPU tasks per submission
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="gpu_tasks_max"
                         :lazy="true"
                         :data="tasksMaxSteps"
                         :disabled="isEmbedding"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ gpu_tasks_max }}
              </span>
            </div>
          </div>

        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 w-full flex-1 overflow-visible pt-0 md:overflow-hidden md:pt-12">

      <div v-if="esmError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="esmError === 'worker'">
          The embedding worker could not be started. Try reloading the page.
        </template>
        <template v-else-if="esmError === 'file_count'">
          Please upload exactly one protein FASTA file.
        </template>
        <template v-else-if="esmError === 'fasta'">
          Not recognised as a protein FASTA. Accepted extensions are
          <b>.faa</b>, <b>.fa</b>, <b>.fasta</b>, <b>.fas</b> and <b>.pep</b>, optionally gzipped.
        </template>
        <template v-else-if="esmError === 'memory'">
          Ran out of memory. Try a smaller file, or split it into several.
        </template>
        <template v-else-if="esmError === 'gpu_lost'">
          The GPU driver dropped the device, and the run could not be retried automatically.
          Try again with <b>GPU acceleration</b> switched off, or with a smaller file: long
          proteins are what drive GPU memory use, since attention cost grows with the square
          of length.
          <span v-if="esmGpuEvent" class="block mt-1 font-mono text-xs">{{ esmGpuEvent }}</span>
        </template>
        <template v-else-if="esmError === 'wasm_panic'">
          The embedding engine stopped unexpectedly. The tab has recovered, so you can try
          again; running without <b>GPU acceleration</b> makes it less likely to recur.
        </template>
        <template v-else>
          {{ esmError }}
        </template>
      </div>

      <div v-if="tabName=='ProteinEmbeddings'">

        <!-- Model loading -->
        <div v-if="isLoadingEsmModel" class="mx-6 mb-4 p-3 bg-blue-50 border border-blue-200 rounded-md">
          <div class="flex items-center gap-2 text-sm">
            <Loader2 class="w-4 h-4 text-blue-500 animate-spin" />
            <span class="font-semibold text-gray-800">
              {{ modelLoadingLabel }}
            </span>
          </div>
        </div>

        <!-- Which backend actually ran -->
        <div v-else-if="esmBackend" class="mx-6 mb-4 text-sm text-gray-500 truncate">
          <span v-if="esmBackend === 'webgpu'">
            GPU<span v-if="activeAdapterName"> · {{ activeAdapterName }}</span>
          </span>
          <template v-else>
            <span>CPU<span v-if="esmBackendFallbackReason"> · GPU unavailable</span></span>
            <span> · considerably slower; consider a smaller file</span>
          </template>
        </div>

        <!-- Upload dropbox - always visible when not embedding -->
        <div v-if="!isEmbedding"
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
            Drop or click to upload your <b>protein FASTA file</b>
          </p>
        </div>

        <div v-else class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
          <Loader2 class="w-6 h-6 text-amber-500 animate-spin"/>
          <p class="text-sm text-gray-500">
            Embedding {{ embeddingFilesArray.join(', ') }}...
          </p>
          <p v-if="embeddingProgress.total > 0" class="text-xs text-gray-400 mt-1">
            Processed proteins: {{ embeddingProgress.done }}/{{ embeddingProgress.total }}
            — {{ progressPercent }}%
          </p>
        </div>

        <!-- Show uploaded files with per-file status -->
        <div v-if="uploadedFileNames.length > 0" class="mx-6 mr-0 mt-4 max-h-48 overflow-y-auto">
          <div v-for="fileName in uploadedFileNames" :key="fileName"
               class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="isFileEmbedded(fileName)" class="w-4 h-4 text-green-500"/>
            <Loader2 v-else-if="isFileEmbedding(fileName)" class="w-4 h-4 text-amber-500 animate-spin"/>
            <span class="flex-grow text-sm font-mono truncate">
              {{ fileName }}
            </span>
          </div>
        </div>

        <!-- Results actions -->
        <div v-if="uploadedFileNames.length > 0" class="mx-6 mr-0 mt-4 flex items-center gap-2">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
          <Button v-if="embeddingsComputed" variant="outline" size="sm" @click="downloadTsv">
            <Download class="mr-2 h-4 w-4" />
            Download TSV
          </Button>
        </div>

        <div v-if="hasProjection" class="px-6 mt-4">
          <EmbeddingScatter
            :coords="esmResult.coords"
            :ids="esmResult.ids"
            :lengths="esmResult.lengths"
          />
        </div>

        <div v-if="embeddingsComputed" class="px-6 mt-4">
          <div class="text-sm text-gray-500 truncate mb-2">
            {{ summaryLine }}
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
import { defineComponent, ref, Ref, computed, onMounted, watch } from "vue";
import { useStore } from "vuex";
import { useDropzone } from "vue3-dropzone";
import { useActions, useState } from "vuex-composition-helpers";
import VueSelect from "vue3-select-component";
import "vue3-select-component/styles";
import VueSlider from "vue-3-slider-component";
import { Check, FileUp, Loader2, Info, ScanBox, Download, Trash2 } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import ProteinEmbeddingsHelpCollapsible from "@/components/help/ProteinEmbeddingsHelpCollapsible.vue";
import DataTable from "@/components/pages/taxonomic-id/DataTable.vue";
import EmbeddingScatter from "@/components/pages/protein-embeddings/EmbeddingScatter.vue";
import { columns, ProteinEmbeddingRow } from "@/components/pages/protein-embeddings/columns";
import { GpuAdapterInfo, ProteinEmbeddingResult } from "@/types";
import { proteinFastaExtensionsWithDotAndCompressList, regExpForAnyProteinFasta } from "@/utils";
import { saveTextFile } from "@/platform/files";
import { isWebGpuAvailable, getWebGLRendererLabel } from "@/platform/gpu";

export default defineComponent({
  name: "ProteinEmbeddingsPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    VueSelect,
    VueSlider,
    Check,
    FileUp,
    Loader2,
    Info,
    ScanBox,
    Download,
    Trash2,
    Button,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    ProteinEmbeddingsHelpCollapsible,
    DataTable,
    EmbeddingScatter
  },
  setup() {
    const store = useStore();
    const uploadedFileNames: Ref<string[]> = ref([]);

    // snake_case for the values handed to the worker; camelCase for UI-only state.
    const use_gpu: Ref<boolean> = ref(false);
    const gpu_power_pref: Ref<number> = ref(1);
    const gpu_tasks_max: Ref<number> = ref(32);
    // Doubling, not linear: the useful range spans two orders of magnitude. 256 is about one
    // forward pass per submission — flat there on a slow GPU, but not on a fast one, where
    // there is less real work for the per-submission cost to hide behind.
    const tasksMaxSteps = [1, 2, 4, 8, 16, 32, 64, 128, 256];
    const auto_device: Ref<boolean> = ref(true);
    const availableAdapters: Ref<GpuAdapterInfo[]> = ref([]);
    const gpuQueryDone: Ref<boolean> = ref(false);
    const gpuSupported: Ref<boolean> = ref(false);

    const { embedFile, resetAllResults_esm, initEmbedderWorkers, listGpuAdapters } = useActions(["embedFile", "resetAllResults_esm", "initEmbedderWorkers", "listGpuAdapters"]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const { allResults_esm } = useState(["allResults_esm"]) as any;

    // The embedder worker is spawned lazily from this page rather than eagerly in App.vue
    // (as the other tabs do), because its wasm module is an order of magnitude larger than
    // the other bridges and would otherwise be loaded for every visitor of every tab.
    onMounted(async () => {
      gpuSupported.value = isWebGpuAvailable();
      if (!store.state.workerState.workers_esm.length) {
        initEmbedderWorkers(1);
      }
      if (gpuSupported.value) {
        // The page re-mounts on every tab switch; the adapter list does not change, so
        // reuse whatever a previous mount already found.
        if (!store.state.gpuAdapters.length) await listGpuAdapters();
        availableAdapters.value = store.state.gpuAdapters;
        gpuQueryDone.value = true;
        gpu_power_pref.value = availableAdapters.value[0]?.index ?? 0;
      }
    });

    watch(auto_device, async (auto) => {
      if (auto || gpuQueryDone.value) return;
      await listGpuAdapters();
      availableAdapters.value = store.state.gpuAdapters;
      gpuQueryDone.value = true;
    });

    // Only meaningful with more than one adapter: with one, its name already is the renderer.
    const rendererHint = computed<string>(() =>
      availableAdapters.value.length > 1 && availableAdapters.value.some(a => !a.identified)
        ? getWebGLRendererLabel()
        : ""
    );

    const activeAdapterName = computed<string>(() =>
      availableAdapters.value.find(a => a.index === gpu_power_pref.value)?.name ?? ""
    );

    /** Options for the adapter VueSelect, including the "nothing to choose" placeholders. */
    const adapterOptions = computed<{ label: string; value: number }[]>(() => {
      if (!gpuSupported.value) return [{ label: "WebGPU not available in this browser", value: 0 }];
      if (!gpuQueryDone.value) return [{ label: "Detecting adapters…", value: 0 }];
      if (!availableAdapters.value.length) {
        return [{ label: "No WebGPU adapter found, CPU will be used", value: 0 }];
      }
      return availableAdapters.value.map(a => ({ label: a.name, value: a.index }));
    });

    // cubecl initialises the device once per wasm instance, so a new value needs a fresh
    // worker; statics reset with the instance, which is all the guard checks.
    watch(gpu_tasks_max, () => {
      if (use_gpu.value && allResults_esm.value.backend === "webgpu") initEmbedderWorkers(1);
    });

    function onDropSample(acceptFiles: File[]): void {
      uploadedFileNames.value = acceptFiles.map(f => f.name);
      embedFile({
        acceptFiles,
        use_gpu: use_gpu.value,
        // Automatic means "not integrated", which is cubecl's DefaultDevice.
        gpu_power_pref: auto_device.value ? 1 : gpu_power_pref.value,
        gpu_tasks_max: gpu_tasks_max.value,
      });
    }

    function resetAll(): void {
      uploadedFileNames.value = [];
      resetAllResults_esm();
    }

    function getSampleNameFromFile(fileName: string): string {
      return fileName.replace(regExpForAnyProteinFasta, "");
    }

    const esmResult = computed<ProteinEmbeddingResult | null>(
      () => allResults_esm.value.result as ProteinEmbeddingResult | null
    );

    const hasProjection = computed<boolean>(() => (esmResult.value?.coords?.length ?? 0) > 0);

    // Every protein: `DataTable` does not virtualise, so a large proteome is slow to render.
    const tableData = computed<ProteinEmbeddingRow[]>(() => {
      const r = esmResult.value;
      if (!r) return [];
      const rows: ProteinEmbeddingRow[] = [];
      for (let i = 0; i < r.nSequences; i++) {
        rows.push({
          id: r.ids[i],
          length: r.lengths[i],
          truncated: r.truncated[i],
          coords: `${r.coords[i * 2].toFixed(2)}, ${r.coords[i * 2 + 1].toFixed(2)}`,
        });
      }
      return rows;
    });


    const embeddingDim = computed<number>(() => esmResult.value?.dim ?? 320);

    /** One-line run summary, in the same idiom AMR detection uses above its table. */
    const summaryLine = computed<string>(() => {
      const r = esmResult.value;
      if (!r) return "";
      const parts = [
        `${r.nSequences} proteins`,
        `${r.dim} dimensions`,
        r.backend === "webgpu" ? "GPU" : "CPU",
        `${(r.elapsedMs / 1000).toFixed(1)} s`,
        r.batchMin === r.batchMax
          ? `batches of ${r.batchMax}`
          : `batches of ${r.batchMin}–${r.batchMax}`,
      ];
      const nTrunc = r.truncated.filter(Boolean).length;
      if (nTrunc > 0) parts.push(`${nTrunc} truncated`);
      return parts.join(" · ");
    });

    /**
     * One row per protein: identifier, length, truncation flag, then `dim` values.
     * toPrecision(7) covers float32's ~7 significant digits, so the round-trip is lossless.
     */
    function embeddingTsv(r: ProteinEmbeddingResult): string {
      const header = ["protein_id", "length", "truncated", "umap_1", "umap_2",
                      ...Array.from({ length: r.dim }, (_, i) => `dim_${i}`)].join("\t");
      const lines = [header];
      for (let i = 0; i < r.nSequences; i++) {
        const row = new Array<string>(r.dim + 5);
        row[0] = r.ids[i];
        row[1] = String(r.lengths[i]);
        row[2] = r.truncated[i] ? "true" : "false";
        // Ahead of the 320 dimensions, so they are visible without scrolling past them.
        row[3] = r.coords[i * 2].toPrecision(7);
        row[4] = r.coords[i * 2 + 1].toPrecision(7);
        const off = i * r.dim;
        for (let j = 0; j < r.dim; j++) {
          row[j + 5] = r.vectors[off + j].toPrecision(7);
        }
        lines.push(row.join("\t"));
      }
      return lines.join("\n") + "\n";
    }

    async function downloadTsv(): Promise<void> {
      const r = esmResult.value;
      if (!r) return;
      await saveTextFile(
        embeddingTsv(r),
        `${r.sampleName}_esm2_embeddings.tsv`,
        "text/tab-separated-values;charset=utf-8",
      );
    }

    const tableColumns = columns;

    const {
      getRootProps: getRootPropsSample,
      getInputProps: getInputPropsSample,
      isDragActive: isDragActiveSample,
      ...restSample
    } = useDropzone({
      onDrop: onDropSample,
      accept: proteinFastaExtensionsWithDotAndCompressList,
      multiple: false
    });

    return {
      use_gpu,
      gpu_power_pref,
      gpu_tasks_max,
      tasksMaxSteps,
      auto_device,
      rendererHint,
      availableAdapters,
      gpuQueryDone,
      gpuSupported,
      activeAdapterName,
      adapterOptions,
      uploadedFileNames,
      resetAll,
      getRootPropsSample,
      getInputPropsSample,
      isDragActiveSample,
      onDropSample,
      getSampleNameFromFile,
      tableColumns,
      tableData,
      embeddingDim,
      summaryLine,
      esmResult,
      hasProjection,
      downloadTsv,
      allResults_esm,
      store,
      ...restSample,
    };
  },
  computed: {
    esmError(): string | null {
      return this.store.getters.esmError;
    },
    isEmbedding(): boolean {
      return this.store.getters.isEmbedding;
    },
    isLoadingEsmModel(): boolean {
      return this.store.getters.isLoadingEsmModel;
    },
    esmModelLoadStage(): string {
      return this.store.getters.esmModelLoadStage;
    },
    modelLoadingLabel(): string {
      if (this.esmModelLoadStage === 'engine') {
        return 'Loading the embedding engine (about 11 MB, first use only)…';
      }
      return this.esmModelLoadStage === 'download'
        ? 'Downloading the ESM-2 model (about 14 MB, first use only)…'
        : `Initialising the model on the ${this.use_gpu ? 'GPU' : 'CPU'}…`;
    },
    esmBackend(): string | null {
      return this.store.getters.esmBackend;
    },
    esmBackendFallbackReason(): string | null {
      return this.store.getters.esmBackendFallbackReason;
    },
    esmGpuEvent(): string | null {
      return this.store.getters.esmGpuEvent;
    },
    embeddingsComputed(): boolean {
      return this.store.getters.embeddingsComputed;
    },
    embeddingFilesSet(): Set<string> {
      return this.store.getters.embeddingFiles;
    },
    embeddingFilesArray(): string[] {
      return Array.from(this.embeddingFilesSet);
    },
    embeddingProgress(): { done: number, total: number } {
      return this.store.getters.embeddingProgress;
    },
    progressPercent(): number {
      const p = this.embeddingProgress;
      return p.total > 0 ? Math.round((p.done / p.total) * 100) : 0;
    },
  },

  methods: {
    clear(): void {
      this.resetAll();
    },
    isFileEmbedded(fileName: string): boolean {
      const sampleName = this.getSampleNameFromFile(fileName);
      return this.allResults_esm.result?.sampleName === sampleName;
    },
    isFileEmbedding(fileName: string): boolean {
      const sampleName = this.getSampleNameFromFile(fileName);
      return this.embeddingFilesSet.has(sampleName);
    }
  },
});
</script>

<style scoped>
</style>
