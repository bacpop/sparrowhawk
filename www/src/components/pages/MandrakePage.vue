<template>
  <div class="flex flex-col gap-6 md:flex-row md:gap-0">
    <div class="w-full md:w-[350px] md:shrink-0">
      <h1 class="text-2xl font-medium flex items-center gap-2 mb-4">
        <ChartScatter class="w-6 h-6" />
        Clustering
      </h1>

      <MandrakeHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <label for="mandrake-labels" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Optional unheadered sample-name/tab/label rows. Every input sample must appear exactly once.</p>
                </TooltipContent>
              </Tooltip>
              Optional sample labels
            </label>
            <input
              id="mandrake-labels"
              ref="labelsInput"
              class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1"
              type="file"
              accept=".tsv,.txt"
              :disabled="isRunning"
              @change="onLabelsChange"
            >
            <p v-if="selectedLabelsFile" class="mt-1 text-xs text-gray-500 font-mono truncate">
              {{ selectedLabelsFile.name }}
            </p>
            <p v-if="labelError" class="mt-1 text-xs text-red-700" role="alert">{{ labelError }}</p>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="mandrake-hdbscan" type="checkbox" v-model="runHdbscan" :disabled="isRunning">
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Run the fixed browser HDBSCAN preset on the final two-dimensional embedding.</p>
              </TooltipContent>
            </Tooltip>
            <label for="mandrake-hdbscan" class="text-sm">Run HDBSCAN after embedding</label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="mandrake-sound" type="checkbox" v-model="soundEnabled">
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Play a tone pair for each embedding frame, pitched by how far the layout moved on each axis.</p>
              </TooltipContent>
            </Tooltip>
            <label for="mandrake-sound" class="text-sm">Add progress sounds</label>
          </div>

          <div v-if="source === 'sketch'" class="flex flex-col gap-4">
            <div>
              <label for="mandrake-sketch-distance" class="flex items-center gap-1 text-sm">
                <Tooltip>
                  <TooltipTrigger as-child>
                    <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                  </TooltipTrigger>
                  <TooltipContent>
                    <p class="max-w-xs">Choose core-distance regression or a selected-k Jaccard distance from the sketch database.</p>
                  </TooltipContent>
                </Tooltip>
                Sketch distance
              </label>
              <select id="mandrake-sketch-distance" v-model="sketchDistance" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" :disabled="isRunning">
                <option value="core" :disabled="sketchKmers.length < 2">Core distance (requires at least two k-mers)</option>
                <option value="jaccard">Jaccard distance</option>
              </select>
            </div>

            <div v-if="sketchDistance === 'jaccard'">
              <label for="mandrake-sketch-kmer" class="flex items-center gap-1 text-sm">
                <Tooltip>
                  <TooltipTrigger as-child>
                    <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                  </TooltipTrigger>
                  <TooltipContent>
                    <p class="max-w-xs">Select one k-mer length stored in the uploaded metadata file.</p>
                  </TooltipContent>
                </Tooltip>
                Jaccard k-mer
              </label>
              <select id="mandrake-sketch-kmer" v-model.number="jaccardKmer" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" :disabled="isRunning || !sketchKmers.length">
                <option v-for="kmer in sketchKmers" :key="kmer" :value="kmer">{{ kmer }}</option>
              </select>
            </div>
          </div>

          <div>
            <label for="mandrake-sparsification" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Choose k-nearest-neighbour or strict normalized distance threshold sparsification.</p>
                </TooltipContent>
              </Tooltip>
              Sparsification
            </label>
            <select id="mandrake-sparsification" v-model="mode" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" :disabled="isRunning || source === 'sketch'">
              <option value="knn">k-nearest neighbours</option>
              <option value="threshold">Distance threshold</option>
            </select>
          </div>

          <div>
            <label for="mandrake-sparsification-value" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">{{ mode === 'knn' ? 'Number of neighbours retained per sample.' : 'Strict normalized distance threshold in the range (0, 1].' }}</p>
                </TooltipContent>
              </Tooltip>
              {{ mode === "knn" ? "Neighbours per sample" : "Distance threshold" }}
            </label>
            <input
              id="mandrake-sparsification-value"
              v-model.number="sparsificationValue"
              class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1"
              type="number"
              :min="mode === 'knn' ? 0 : 0.001"
              :max="mode === 'knn' ? 10000 : 1"
              :step="mode === 'knn' ? 1 : 0.01"
              :disabled="isRunning"
            >
          </div>

          <div>
            <label for="mandrake-perplexity" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Conditional-probability perplexity (local vs global attention), in the inclusive range 5 to 100.</p>
                </TooltipContent>
              </Tooltip>
              Perplexity
            </label>
            <input id="mandrake-perplexity" v-model.number="perplexity" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" type="number" min="5" max="100" step="1" :disabled="isRunning">
          </div>

          <div>
            <label for="mandrake-max-updates" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Target number of stochastic update attempts (increase for larger numbers of samples).</p>
                </TooltipContent>
              </Tooltip>
              Maximum updates
            </label>
            <input id="mandrake-max-updates" v-model.number="maxUpdates" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" type="number" min="1" step="1000000" :disabled="isRunning">
          </div>

          <div>
            <label for="mandrake-repulsion" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Repulsion samples used per update attempt.</p>
                </TooltipContent>
              </Tooltip>
              Repulsion samples
            </label>
            <input id="mandrake-repulsion" v-model.number="repulsionSamples" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" type="number" min="1" step="1" :disabled="isRunning">
          </div>

          <div>
            <label for="mandrake-learning-rate" class="flex items-center gap-1 text-sm">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Initial learning rate for the stochastic embedding (rate at which embedding changes).</p>
                </TooltipContent>
              </Tooltip>
              Learning rate
            </label>
            <input id="mandrake-learning-rate" v-model.number="learningRate" class="w-full border border-gray-300 rounded-md text-sm p-2 mt-1" type="number" min="0.0001" step="0.1" :disabled="isRunning">
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="mandrake-exaggeration" type="checkbox" v-model="initialExaggeration" :disabled="isRunning">
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Apply initial attraction exaggeration, which quadruples attraction forces in the first 10%.</p>
              </TooltipContent>
            </Tooltip>
            <label for="mandrake-exaggeration" class="text-sm">Use initial exaggeration</label>
          </div>
        </div>
      </TooltipProvider>

      <div class="flex flex-wrap gap-2 mt-6">
        <Button :disabled="!canRun || isRunning" @click="runEmbedding">
          {{ isRunning ? "Embedding…" : "Run Mandrake" }}
        </Button>
        <Button v-if="isRunning" variant="outline" @click="cancel">Cancel</Button>
      </div>
    </div>

    <div class="min-w-0 w-full flex-1 overflow-visible pt-0 md:overflow-hidden md:pt-12">
      <div v-if="errorMessage" class="mx-6 mr-0 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800" role="alert">
        {{ errorMessage }}
      </div>

      <div
        v-bind="getRootProps()"
        :class="[
          'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
          isRunning ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer hover:border-gray-400'
        ]"
      >
        <input v-bind="getInputProps()" :disabled="isRunning" />
        <FileUp class="w-6 h-6" />
        <p v-if="isDragActive">Drop files here …</p>
        <p v-else>Drop or click to upload an <b>alignment, accessory table, or sketch database</b></p>
        <p class="text-xs text-gray-400">FASTA, Rtab/TSV, or paired .skm/.skd files</p>
      </div>

      <div v-if="selectedFile || source === 'sketch'" class="mx-6 mr-0 mt-4">
        <div v-if="selectedFile" class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
          <Check v-if="!isRunning" class="w-4 h-4 text-green-500" />
          <Loader2 v-else class="w-4 h-4 text-blue-500 animate-spin" />
          <span class="flex-grow text-sm font-mono truncate">{{ selectedFile.name }}</span>
          <span class="text-xs text-gray-500">{{ sourceLabel }}</span>
        </div>
        <template v-if="source === 'sketch'">
          <div class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="sketchMetadataFile && !sketchMetadataLoading" class="w-4 h-4 text-green-500" />
            <Loader2 v-else-if="sketchMetadataLoading" class="w-4 h-4 text-blue-500 animate-spin" />
            <span class="flex-grow text-sm font-mono truncate">{{ sketchMetadataFile?.name ?? "No .skm metadata file" }}</span>
            <span class="text-xs text-gray-500">Sketch metadata</span>
          </div>
          <div class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="sketchDataFile" class="w-4 h-4 text-green-500" />
            <span class="flex-grow text-sm font-mono truncate">{{ sketchDataFile?.name ?? "No .skd data file" }}</span>
            <span class="text-xs text-gray-500">Sketch data</span>
          </div>
          <p v-if="sketchMetadataLoading" class="text-xs text-gray-500">Reading sketch metadata …</p>
        </template>
        <p v-if="inputError" class="mt-1 text-xs text-red-700" role="alert">{{ inputError }}</p>
        <Button variant="outline" size="sm" class="mt-2" :disabled="isRunning" @click="clearAll">
          <Trash2 class="mr-1 h-3 w-3" />
          Clear input
        </Button>
      </div>

      <div v-else-if="inputError" class="mx-6 mr-0 mt-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800" role="alert">
        {{ inputError }}
      </div>

      <div v-if="sketchLoading" class="mx-6 mr-0 mt-4 p-3 bg-blue-50 border border-blue-200 rounded-md" role="status">
        <div class="flex items-center gap-2 text-sm">
          <Loader2 class="w-4 h-4 text-blue-500 animate-spin" />
          <span class="font-semibold text-gray-800">Sketch distance phase</span>
        </div>
        <p class="mt-1 ml-6 text-xs text-gray-600">Loading the paired sketch database and calculating distances</p>
      </div>

      <div v-if="isRunning || distanceProgress.maximum > 0" class="mx-6 mr-0 mt-4 p-3 bg-blue-50 border border-blue-200 rounded-md">
        <div class="flex items-center justify-between gap-2 text-sm">
          <span class="font-semibold text-gray-800">Distance phase</span>
          <span>{{ distancePercent }}%</span>
        </div>
        <div class="h-2 mt-2 overflow-hidden rounded-full bg-blue-100" aria-hidden="true">
          <div class="h-full rounded-full bg-blue-500 transition-[width] duration-150" :style="{ width: `${distancePercent}%` }" />
        </div>
        <p class="mt-2 text-xs text-gray-600 font-mono">
          {{ distanceProgress.completed.toLocaleString() }} / {{ distanceProgress.maximum.toLocaleString() }} rows
        </p>
      </div>

      <div v-if="isRunning || embeddingProgress.maximum > 0" class="mx-6 mr-0 mt-4 p-3 bg-blue-50 border border-blue-200 rounded-md">
        <div class="flex items-center justify-between gap-2 text-sm">
          <span class="font-semibold text-gray-800">{{ isRunning ? "Embedding phase" : "Embedding complete" }}</span>
          <span>{{ embeddingPercent }}%</span>
        </div>
        <div class="h-2 mt-2 overflow-hidden rounded-full bg-blue-100" aria-hidden="true">
          <div class="h-full rounded-full bg-orange-500 transition-[width] duration-150" :style="{ width: `${embeddingPercent}%` }" />
        </div>
        <p class="mt-2 text-xs text-gray-600 font-mono">
          {{ embeddingProgress.completed.toLocaleString() }} / {{ embeddingProgress.maximum.toLocaleString() }} updates
          <span v-if="Number.isFinite(embeddingProgress.eq)"> · Eq {{ embeddingProgress.eq.toFixed(4) }}</span>
        </p>
      </div>

      <div v-if="clustering" class="mx-6 mr-0 mt-4 p-3 bg-amber-50 border border-amber-300 rounded-md" role="status">
        <div class="flex items-center gap-2 text-sm">
          <Loader2 class="w-4 h-4 text-amber-500 animate-spin" />
          <span class="font-semibold text-gray-800">HDBSCAN labelling</span>
        </div>
        <p class="mt-1 ml-6 text-xs text-gray-600">Labelling the final embedding</p>
      </div>

      <div v-if="liveEmbedding" class="mx-6 mr-0 mt-4 border border-gray-200 rounded-md bg-white p-4">
        <div class="flex flex-col gap-3 sm:flex-row sm:items-start sm:justify-between">
          <div>
            <h2 class="text-lg font-medium">
              {{ isRunning ? (clustering ? "Labelling clusters …" : "Embedding in progress") : "Final embedding" }}
            </h2>
            <p class="text-sm text-gray-500 mt-1">
              {{ sampleNames.length.toLocaleString() }} samples · two dimensions
              <span v-if="hdbscanLabels"> · {{ hdbscanSummary }}</span>
            </p>
          </div>
          <div v-if="result" class="flex flex-wrap gap-2">
            <Button variant="outline" size="sm" @click="downloadEmbedding"><Download class="mr-1 h-3 w-3" />Embedding</Button>
            <Button variant="outline" size="sm" @click="downloadNames"><Download class="mr-1 h-3 w-3" />Names</Button>
            <Button v-if="hdbscanLabels" variant="outline" size="sm" @click="downloadClusters"><Download class="mr-1 h-3 w-3" />Clusters</Button>
          </div>
        </div>
        <div v-if="clusterError" class="mt-3 p-3 bg-amber-50 border border-amber-300 rounded-md text-sm text-amber-800" role="status">
          HDBSCAN could not label this embedding: {{ clusterError }}
        </div>
        <div v-if="hasBothLabels" class="flex flex-wrap items-center gap-2 mt-4 text-sm" role="group" aria-label="Plot colours">
          <span class="font-medium text-gray-600">Colour by</span>
          <Button variant="outline" size="sm" :class="colourMode === 'manual' ? 'bg-gray-900 text-white hover:bg-gray-800' : ''" @click="colourMode = 'manual'">Manual labels</Button>
          <Button variant="outline" size="sm" :class="colourMode === 'clusters' ? 'bg-gray-900 text-white hover:bg-gray-800' : ''" @click="colourMode = 'clusters'">HDBSCAN clusters</Button>
        </div>
        <div class="mt-4 min-h-[440px] overflow-hidden rounded-md border border-gray-200 bg-[#fcfcfb]">
          <MandrakeEmbeddingPlot
            :embedding="liveEmbedding"
            :names="sampleNames"
            :labels="activeLabels ?? undefined"
            :noise-label="colourMode === 'clusters' ? 'Noise' : undefined"
            :run-key="runKey"
          />
        </div>
      </div>

      <div v-else-if="!isRunning && !errorMessage" class="mx-6 mr-0 mt-4 min-h-[420px] border border-gray-200 rounded-md bg-gray-50 flex flex-col justify-center items-center gap-3 p-8 text-center">
        <ChartScatter class="w-10 h-10 text-gray-400" />
        <h2 class="text-lg font-medium text-gray-800">Your embedding will appear here</h2>
        <p class="text-sm text-gray-500">Choose an input file and parameters, then run Mandrake.</p>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import { computed, onBeforeUnmount, ref, shallowRef, watch } from "vue";
import { useDropzone } from "vue3-dropzone";
import { ChartScatter, Check, Download, FileUp, Info, Loader2, Trash2 } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import MandrakeHelpCollapsible from "@/components/help/MandrakeHelpCollapsible.vue";
import MandrakeEmbeddingPlot from "@/components/MandrakeEmbeddingPlot.vue";
import { saveTextFile } from "@/platform/files";
import { MandrakeSonifier } from "@/lib/mandrakeSound";
import {
  MandrakeRunner,
  type MandrakeDistanceProgress,
  type MandrakeProgress,
  type MandrakeResult,
  type MandrakeSketchFiles,
  type MandrakeSettings,
  type MandrakeUpdate,
} from "@/workers/Mandrake";

type InputSource = "alignment" | "accessory" | "sketch";

defineProps<{ tabName: string }>();

const mandrakeInputExtensions = [
  ".fa", ".fasta", ".fas", ".fna",
  ".fa.gz", ".fasta.gz", ".fas.gz", ".fna.gz",
  ".rtab", ".Rtab", ".tsv",
  ".rtab.gz", ".Rtab.gz", ".tsv.gz",
  ".skm", ".skd", ".gz",
];

const labelsInput = ref<HTMLInputElement | null>(null);
const source = ref<InputSource | null>(null);
const selectedFile = ref<File | null>(null);
const sketchMetadataFile = ref<File | null>(null);
const sketchDataFile = ref<File | null>(null);
const sketchKmers = ref<number[]>([]);
const sketchMetadataLoading = ref(false);
const sketchDistance = ref<"core" | "jaccard">("core");
const jaccardKmer = ref<number | null>(null);
const selectedLabelsFile = ref<File | null>(null);
const inputError = ref("");
const labelError = ref("");
const mode = ref<"knn" | "threshold">("knn");
const sparsificationValue = ref(15);
const perplexity = ref(30);
const maxUpdates = ref(1_000_000);
const repulsionSamples = ref(5);
const learningRate = ref(1);
const initialExaggeration = ref(false);
const runHdbscan = ref(false);
const soundEnabled = ref(false);
const sonifier = new MandrakeSonifier();
const isRunning = ref(false);
const errorMessage = ref("");
const result = ref<MandrakeResult | null>(null);
const liveEmbedding = shallowRef<Float64Array | null>(null);
const sampleNames = ref<string[]>([]);
const labels = ref<string[] | null>(null);
const hdbscanLabels = shallowRef<Int32Array | null>(null);
const colourMode = ref<"manual" | "clusters">("manual");
const clustering = ref(false);
const sketchLoading = ref(false);
const clusterError = ref("");
const labelContents = ref<string | null>(null);
const distanceProgress = ref<MandrakeDistanceProgress>({ completed: 0, maximum: 0, complete: false });
const embeddingProgress = ref<MandrakeProgress>({ completed: 0, maximum: 0, eq: Number.NaN, complete: false });
const runKey = ref(0);
const runner = new MandrakeRunner();

const { getRootProps, getInputProps, isDragActive } = useDropzone({
  onDrop: chooseFiles,
  accept: mandrakeInputExtensions.join(","),
  multiple: true,
});

const sourceLabel = computed(() => source.value === "alignment"
  ? "Alignment"
  : source.value === "accessory"
    ? "Accessory table"
    : "Sketch database");

const canRun = computed(() => source.value === "sketch"
  ? sketchMetadataFile.value !== null && sketchDataFile.value !== null
  : selectedFile.value !== null && source.value !== null);

const distancePercent = computed(() => {
  if (!distanceProgress.value.maximum) return 0;
  return Math.min(100, Math.round((distanceProgress.value.completed / distanceProgress.value.maximum) * 100));
});

const embeddingPercent = computed(() => {
  if (!embeddingProgress.value.maximum) return 0;
  return Math.min(100, Math.round((embeddingProgress.value.completed / embeddingProgress.value.maximum) * 100));
});

const hdbscanPlotLabels = computed<string[] | null>(() => {
  const clusterLabels = hdbscanLabels.value;
  if (!clusterLabels || clusterLabels.length !== sampleNames.value.length) return null;
  return Array.from(clusterLabels, (label) => label < 0 ? "Noise" : `Cluster ${label}`);
});

const hasBothLabels = computed(() => labels.value !== null && hdbscanPlotLabels.value !== null);
const activeLabels = computed(() => colourMode.value === "clusters" ? hdbscanPlotLabels.value : labels.value);

const clusterCount = computed(() => {
  const clusterLabels = hdbscanLabels.value;
  if (!clusterLabels) return 0;
  return new Set(Array.from(clusterLabels).filter((label) => label >= 0)).size;
});

const hdbscanSummary = computed(() => clusterCount.value === 0
  ? "No HDBSCAN clusters found"
  : `${clusterCount.value} HDBSCAN cluster${clusterCount.value === 1 ? "" : "s"}`);

watch(mode, (nextMode) => {
  sparsificationValue.value = nextMode === "knn" ? 15 : 0.5;
});

watch(source, (nextSource) => {
  if (nextSource === "sketch") mode.value = "knn";
});

watch(soundEnabled, (enabled) => {
  if (enabled) sonifier.enable();
  else sonifier.disable();
});

function detectSource(filename: string): Exclude<InputSource, "sketch"> | null {
  const lowerFilename = filename.toLowerCase();
  const sourceFilename = lowerFilename.endsWith(".gz") ? lowerFilename.slice(0, -3) : lowerFilename;
  const suffix = sourceFilename.slice(sourceFilename.lastIndexOf("."));
  if ([".fa", ".fasta", ".fas", ".fna"].includes(suffix)) return "alignment";
  if ([".rtab", ".tsv"].includes(suffix)) return "accessory";
  return null;
}

function detectSketchKind(filename: string): "metadata" | "data" | null {
  const lowerFilename = filename.toLowerCase();
  if (lowerFilename.endsWith(".skm")) return "metadata";
  if (lowerFilename.endsWith(".skd")) return "data";
  return null;
}

function resetResultState(): void {
  selectedFile.value = null;
  result.value = null;
  liveEmbedding.value = null;
  sampleNames.value = [];
  labels.value = null;
  hdbscanLabels.value = null;
  colourMode.value = "manual";
  clustering.value = false;
  sketchLoading.value = false;
  clusterError.value = "";
  labelContents.value = null;
  distanceProgress.value = { completed: 0, maximum: 0, complete: false };
  embeddingProgress.value = { completed: 0, maximum: 0, eq: Number.NaN, complete: false };
  errorMessage.value = "";
  inputError.value = "";
  labelError.value = "";
}

function clearInputSelection(): void {
  selectedFile.value = null;
  sketchMetadataFile.value = null;
  sketchDataFile.value = null;
  sketchKmers.value = [];
  sketchMetadataLoading.value = false;
  sketchDistance.value = "core";
  jaccardKmer.value = null;
  source.value = null;
}

function clearAll(): void {
  runner.cancel();
  isRunning.value = false;
  clearInputSelection();
  selectedLabelsFile.value = null;
  resetResultState();
}

async function inspectSketchMetadata(file: File): Promise<void> {
  sketchMetadataLoading.value = true;
  sketchKmers.value = [];
  jaccardKmer.value = null;
  try {
    const kmers = await runner.inspectSketchKmers(file);
    if (sketchMetadataFile.value !== file) return;
    sketchKmers.value = kmers;
    jaccardKmer.value = kmers[0] ?? null;
    sketchDistance.value = kmers.length >= 2 ? "core" : "jaccard";
  } catch (error) {
    if (sketchMetadataFile.value === file) inputError.value = error instanceof Error ? error.message : String(error);
  } finally {
    if (sketchMetadataFile.value === file) sketchMetadataLoading.value = false;
  }
}

function chooseFiles(files: File[]): void {
  if (!files.length || isRunning.value) return;
  const sketchKinds = files.map((file) => detectSketchKind(file.name));
  const hasSketch = sketchKinds.some((kind) => kind !== null);
  if (hasSketch) {
    if (files.some((file, index) => sketchKinds[index] === null)) {
      clearInputSelection();
      resetResultState();
      inputError.value = "Sketch input requires one .skm metadata file and one .skd data file; do not mix sketch and sequence inputs.";
      return;
    }
    if (source.value !== "sketch") {
      resetResultState();
      clearInputSelection();
      source.value = "sketch";
    } else {
      resetResultState();
    }
    const metadataFiles = files.filter((_, index) => sketchKinds[index] === "metadata");
    const dataFiles = files.filter((_, index) => sketchKinds[index] === "data");
    if (metadataFiles.length > 1 || dataFiles.length > 1) {
      inputError.value = "Select at most one .skm metadata file and one .skd data file.";
      return;
    }
    const metadataFile = metadataFiles[0];
    const dataFile = dataFiles[0];
    if (metadataFile && sketchMetadataFile.value) {
      inputError.value = "Only one .skm metadata file can be selected.";
    } else if (dataFile && sketchDataFile.value) {
      inputError.value = "Only one .skd data file can be selected.";
    } else {
      inputError.value = "";
      if (metadataFile) {
        sketchMetadataFile.value = metadataFile;
        void inspectSketchMetadata(metadataFile);
      }
      if (dataFile) sketchDataFile.value = dataFile;
    }
    return;
  }

  if (files.length !== 1) {
    clearInputSelection();
    resetResultState();
    inputError.value = "Select one FASTA or Rtab/TSV file, or a paired .skm/.skd database.";
    return;
  }
  const file = files[0];
  const detectedSource = detectSource(file.name);
  clearInputSelection();
  resetResultState();
  if (!detectedSource) {
    inputError.value = "Unsupported input suffix. Use FASTA or Rtab/TSV, optionally followed by .gz, or select one .skm and one .skd file.";
    return;
  }
  selectedFile.value = file;
  source.value = detectedSource;
}

function onLabelsChange(event: Event): void {
  const input = event.target as HTMLInputElement;
  selectedLabelsFile.value = input.files?.[0] ?? null;
  labels.value = null;
  labelContents.value = null;
  labelError.value = "";
  input.value = "";
}

function parseLabels(contents: string, names: string[]): string[] {
  const lines = contents.replace(/\r\n/g, "\n").split("\n");
  if (lines.at(-1) === "") lines.pop();
  const labelsByName = new Map<string, string>();
  lines.forEach((line, index) => {
    const fields = line.split("\t");
    if (fields.length !== 2) throw new Error(`label file row ${index + 1} must contain exactly two tab-separated fields`);
    const [name, label] = fields;
    if (!name) throw new Error(`label file row ${index + 1} has an empty sample name`);
    if (labelsByName.has(name)) throw new Error(`label file contains duplicate sample name: ${name}`);
    labelsByName.set(name, label);
  });
  const nameSet = new Set(names);
  const missing = names.filter((name) => !labelsByName.has(name));
  const extra = Array.from(labelsByName.keys()).filter((name) => !nameSet.has(name));
  if (missing.length || extra.length) {
    const details = [];
    if (missing.length) details.push(`missing names: ${missing.join(", ")}`);
    if (extra.length) details.push(`unknown names: ${extra.join(", ")}`);
    throw new Error(`label/name mismatch (${details.join("; ")})`);
  }
  return names.map((name) => labelsByName.get(name)!);
}

function handleUpdate(update: MandrakeUpdate): void {
  if (update.phase === "sketch") {
    sketchLoading.value = true;
    return;
  }
  if (update.phase === "distance") {
    sketchLoading.value = false;
    distanceProgress.value = update.progress;
    if (update.names?.length) {
      sampleNames.value = update.names;
      if (labelContents.value !== null && labels.value === null) {
        try {
          labels.value = parseLabels(labelContents.value, update.names);
        } catch (error) {
          labelError.value = error instanceof Error ? error.message : String(error);
          errorMessage.value = labelError.value;
          runner.cancel();
        }
      }
    }
    return;
  }
  if (update.phase === "embedding") {
    embeddingProgress.value = update.progress;
    return;
  }
  if (update.phase === "clustering") {
    clustering.value = true;
    return;
  }
  if (soundEnabled.value) sonifier.playFrame(update.embedding);
  liveEmbedding.value = update.embedding;
}

function settings(): MandrakeSettings {
  return {
    mode: source.value === "sketch" ? "knn" : mode.value,
    value: Number(sparsificationValue.value),
    perplexity: Number(perplexity.value),
    maxUpdates: Number(maxUpdates.value),
    repulsionSamples: Number(repulsionSamples.value),
    learningRate: Number(learningRate.value),
    initialExaggeration: initialExaggeration.value,
    hdbscan: runHdbscan.value,
    sketchDistance: sketchDistance.value,
    jaccardKmer: jaccardKmer.value ?? 0,
  };
}

async function runEmbedding(): Promise<void> {
  if (!canRun.value || !source.value) return;
  isRunning.value = true;
  errorMessage.value = "";
  labelError.value = "";
  result.value = null;
  liveEmbedding.value = null;
  sampleNames.value = [];
  labels.value = null;
  hdbscanLabels.value = null;
  colourMode.value = "manual";
  clustering.value = false;
  sketchLoading.value = false;
  clusterError.value = "";
  labelContents.value = null;
  distanceProgress.value = { completed: 0, maximum: 0, complete: false };
  embeddingProgress.value = { completed: 0, maximum: Number(maxUpdates.value), eq: Number.NaN, complete: false };
  runKey.value += 1;
  if (soundEnabled.value) sonifier.enable();
  try {
    labelContents.value = selectedLabelsFile.value ? await selectedLabelsFile.value.text() : null;
    if (source.value === "sketch") {
      if (!sketchMetadataFile.value || !sketchDataFile.value) return;
      const files: MandrakeSketchFiles = { metadata: sketchMetadataFile.value, data: sketchDataFile.value };
      result.value = await runner.runSketch(files, settings(), handleUpdate);
    } else if (selectedFile.value) {
      result.value = await runner.run(source.value, selectedFile.value, settings(), handleUpdate);
    }
    if (!result.value) return;
    sampleNames.value = result.value.names;
    liveEmbedding.value = result.value.embedding;
    hdbscanLabels.value = result.value.hdbscanLabels;
    clusterError.value = result.value.hdbscanError ?? "";
    clustering.value = false;
    colourMode.value = result.value.hdbscanLabels && labels.value === null ? "clusters" : "manual";
  } catch (error) {
    if (!(error instanceof Error) || error.message !== "Mandrake operation cancelled") {
      errorMessage.value = error instanceof Error ? error.message : String(error);
    }
  } finally {
    isRunning.value = false;
  }
}

function cancel(): void {
  runner.cancel();
  isRunning.value = false;
  sketchLoading.value = false;
}

function outputPrefix(): string {
  if (source.value === "sketch") {
    return (sketchMetadataFile.value?.name ?? sketchDataFile.value?.name ?? "mandrake")
      .replace(/\.skm$/i, "")
      .replace(/\.skd$/i, "");
  }
  const filename = selectedFile.value?.name ?? "mandrake";
  return filename.replace(/\.(?:gz)$/i, "").replace(/\.[^.]+$/, "");
}

function downloadEmbedding(): void {
  if (!result.value) return;
  const rows: string[] = [];
  for (let index = 0; index < result.value.embedding.length; index += 2) {
    rows.push(`${result.value.embedding[index].toExponential(17)}\t${result.value.embedding[index + 1].toExponential(17)}`);
  }
  void saveTextFile(`${rows.join("\n")}\n`, `${outputPrefix()}.embedding.txt`);
}

function downloadNames(): void {
  if (!result.value) return;
  void saveTextFile(`${result.value.names.join("\n")}\n`, `${outputPrefix()}.names.txt`);
}

function csvEscape(value: string | number): string {
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function downloadClusters(): void {
  const clusterLabels = hdbscanLabels.value;
  if (!result.value || !clusterLabels || clusterLabels.length !== result.value.names.length) return;
  const rows = ["id,hdbscan_cluster__autocolour"];
  result.value.names.forEach((name, index) => rows.push(`${csvEscape(name)},${clusterLabels[index]}`));
  void saveTextFile(`${rows.join("\n")}\n`, `${outputPrefix()}.embedding_hdbscan_clusters.csv`, "text/csv");
}

onBeforeUnmount(() => {
  runner.cancel();
  sonifier.dispose();
});
</script>
