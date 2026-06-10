<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <Dna class="w-6 h-6" />
        Gene calling
      </h1>

      <GeneCallingHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Number of parallel workers used to process files simultaneously. Higher values speed up batch processing but use more memory.</p>
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
                         :disabled="callingGenes"
                         @change="onNumWorkersChange"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ numWorkers }}
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
                  <p class="max-w-xs">You can select a different translation table if you prefer from NCBI ones.</p>
                </TooltipContent>
              </Tooltip>
              Translation table
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSelect class="flex-grow"
                         v-model="tt"
                         :options="[
                          { label: 'Default/Auto', value: 0},
                          { label: '1. Standard', value: 1},
                          { label: '2. Vertebrate mitochondrial', value: 2},
                          { label: '3. Yeast mitochondrial', value: 3},
                          { label: '4. Mold, protozoan, coelenterate mitochondrial, and mycoplasma/spiroplasma', value: 4},
                          { label: '5. Invertebrate mitochondrial', value: 5},
                          { label: '6. Ciliate, dasycladacean, and hexamita nuclear', value: 6},
                          { label: '9. Echinoderm, and flatworm mitochondrial', value: 9},
                          { label: '10. Euplotid nuclear', value: 10},
                          { label: '11. Bacterial, archaeal, and plant plastid', value: 11},
                          { label: '12. Alternative yeast nuclear', value: 12},
                          { label: '13. Ascidian mitochondrial', value: 13},
                          { label: '14. Alternative flatworm mitochondrial', value: 14},
                          { label: '15. Blepharisma nuclear', value: 15},
                          { label: '16. Chlorophycean mitochondrial', value: 16},
                          { label: '21. Trematode mitochondrial', value: 21},
                          { label: '22. Scenedesmus obliquus mitochondrial', value: 22},
                          { label: '23. Thraustochytrium mitochondrial', value: 23},
                          { label: '24. Rhabdopleuridate mitochondrial', value: 24},
                          { label: '25. Candidate division SR1 and gracilibacteria', value: 25},
                         ]"
                         :isDisabled="callingGenes"
                         :isClearable=false
              />
            </div>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="metag" type="checkbox" v-model="metag" :disabled="callingGenes"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Recommended for FASTA files with lots of short sequences, not necessarily related between them.</p>
              </TooltipContent>
            </Tooltip>
            <label for="metag">
              Use metagenomic mode
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="closed_ends" type="checkbox" v-model="closed_ends" :disabled="callingGenes"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Ignores any gene that might have run off the edge of a contig</p>
              </TooltipContent>
            </Tooltip>
            <label for="closed_ends">
              Ignore truncated genes
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="mask" type="checkbox" v-model="mask" :disabled="callingGenes"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">When finding a gap, or set of unknown (N) nucleotides, we will not bridge over them to call a gene.</p>
              </TooltipContent>
            </Tooltip>
            <label for="mask">
              Break calling on N subsequences
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="non_sd" type="checkbox" v-model="non_sd" :disabled="callingGenes"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Force the algorithm to not use Shine-Dalgarno sequences for calling genes</p>
              </TooltipContent>
            </Tooltip>
            <label for="non_sd">
              Ignore Shine-Dalgarno sequences
            </label>
          </div>

          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Minimum fraction of exact-gene diagnostic k-mers required to annotate a called CDS as an AMR gene.</p>
                </TooltipContent>
              </Tooltip>
              Minimum AMR gene fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_gene_fraction"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.01"
                         :disabled="callingGenes"
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
                  <p class="max-w-xs">Minimum fraction of diagnostic k-mers required to annotate a called CDS as a collapsed AMR gene-group or hierarchy report unit.</p>
                </TooltipContent>
              </Tooltip>
              Minimum AMR gene-group fraction
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="min_gene_group_fraction"
                         :lazy="true"
                         :min="0"
                         :max="1"
                         :interval="0.01"
                         :disabled="callingGenes"
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

      <div v-if="orphosError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="orphosError === 'memory'">
          Error during processing — most likely a memory issue. Try with fewer or smaller files.
        </template>
        <template v-else>
          Error during processing: {{ orphosError }}
        </template>
      </div>

      <div v-if="tabName=='GeneCalling'">

        <!-- Upload dropzone - always visible -->
        <div v-bind='getRootPropsGenome()'
             :class="[
               'p-6 mx-6 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind='getInputPropsGenome()' />
          <FileUp/>
          <p v-if='isDragActiveGenome'>
            Drop files here ...
          </p>
          <p v-else>
            Drop or click to upload your <b>FASTA file(s)</b>
          </p>
        </div>

        <!-- Per-file in-flight status list -->
        <div v-if="uploadedFileNames.length > 0" class="mx-6 mt-4 max-h-48 overflow-y-auto">
          <div v-for="fileName in uploadedFileNames" :key="fileName"
               class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Check v-if="isFileDone(fileName)" class="w-4 h-4 text-green-500"/>
            <Loader2 v-else-if="isFileInFlight(fileName)" class="w-4 h-4 text-amber-500 animate-spin"/>
            <span class="flex-grow text-sm font-mono truncate">{{ fileName }}</span>
          </div>
        </div>

        <!-- Progress block -->
        <div v-if="callingGenes" class="mx-6 mt-4 p-3 bg-blue-50 border border-blue-200 rounded-md">
          <div class="flex items-center gap-2 text-sm">
            <Loader2 class="w-4 h-4 text-blue-500 animate-spin" />
            <span v-if="geneCallingProgressTotal <= 1" class="font-semibold text-gray-800">
              {{ geneCallingStepLabel }}
            </span>
            <span v-else class="text-gray-800">
              Processed files: {{ resultCount }}/{{ geneCallingProgressTotal }}
              — {{ Math.round(resultCount / geneCallingProgressTotal * 100) }}%
            </span>
          </div>
        </div>

        <!-- Action bar -->
        <div v-if="genesCalled || callingGenes" class="mx-6 mt-4 flex gap-2">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
          <template v-if="resultCount >= 2">
            <Button variant="outline" size="sm" @click="downloadZip">
              <Download class="mr-2 h-4 w-4" />
              Download all (.zip)
            </Button>
            <Button variant="outline" size="sm" @click="downloadTarGz">
              <Download class="mr-2 h-4 w-4" />
              Download all (.tar.gz)
            </Button>
          </template>
        </div>

        <!-- Results table -->
        <div v-if="genesCalled" class="mx-6 mt-4 overflow-x-auto">
          <table class="w-full text-sm border border-gray-200 rounded-md">
            <thead>
              <tr class="bg-gray-50 text-left">
                <th class="px-3 py-2 font-medium text-gray-700">File name</th>
                <th class="px-3 py-2 font-medium text-gray-700">Genes called</th>
                <th class="px-3 py-2 font-medium text-gray-700">AMR hits</th>
                <th class="px-3 py-2 font-medium text-gray-700">Download</th>
                <th class="px-3 py-2 font-medium text-gray-700">View</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(result, fileName) in orphosResults" :key="fileName"
                  class="border-t border-gray-100">
                <td class="px-3 py-2 font-mono truncate max-w-[200px]">{{ result.fileName }}</td>
                <td class="px-3 py-2">{{ result.geneCount }}</td>
                <td class="px-3 py-2">
                  <span>{{ result.amrHitCount }}</span>
                  <span v-if="result.amrError" class="ml-2 text-xs text-amber-700">AMR warning</span>
                </td>
                <td class="px-3 py-2">
                  <div class="flex flex-wrap gap-2">
                    <Button variant="outline" size="sm" @click="downloadGff(result)">
                      <Download class="mr-1 h-3 w-3" />
                      .gff
                    </Button>
                    <Button v-if="result.amrAnnotationTsv" variant="outline" size="sm" @click="downloadAmrTsv(result)">
                      <Download class="mr-1 h-3 w-3" />
                      .amr.tsv
                    </Button>
                  </div>
                </td>
                <td class="px-3 py-2">
                  <Button variant="outline" size="sm"
                          :class="selectedGenome === result.fileName ? 'bg-blue-50' : ''"
                          @click="selectGenome(result.fileName)">
                    <Eye class="mr-1 h-3 w-3" />
                    View
                  </Button>
                </td>
              </tr>
            </tbody>
          </table>
        </div>

        <!-- Genome viewer -->
        <div v-if="selectedGenome" class="mx-6 mt-4 overflow-hidden" style="height: 600px;">
          <div ref="genomeViewerContainer" style="height: 600px;"></div>
        </div>

      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, Ref, computed, watch, nextTick, onBeforeUnmount } from "vue";
import { useStore } from "vuex";
import VueSelect from "vue3-select-component";
import "vue3-select-component/styles";
import VueSlider from 'vue-3-slider-component';
import { useDropzone } from "vue3-dropzone";
import { useActions } from "vuex-composition-helpers";
import { Check, FileUp, Loader2, Info, Dna, Download, Eye, Trash2 } from "lucide-vue-next";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import GeneCallingHelpCollapsible from "@/components/help/GeneCallingHelpCollapsible.vue";
import { fastaExtensionsWithDotAndCompressList } from "@/utils";
import { GeneCallResult } from "@/types";
import { buildZip, buildTarGz } from "@/archiveUtils";
import * as ReactDOM from 'react-dom/client';
import * as React from 'react';
import { GeneCallingGenomeViewer } from '@/components/gene-calling/GeneCallingGenomeViewer';

const STEP_LABELS: Record<string, string> = {
    creating_interface:   'Creating interface',
    reading_fasta:        'Reading FASTA file',
    compressing_indexing: 'Compressing and indexing',
    calling_genes:        'Calling genes',
    genes_called:         'Genes called',
    writing_gff:          'Writing and indexing GFF file',
    detecting_amr:        'Detecting AMR in called genes',
    writing_annotated_gff:'Writing annotated GFF file',
    done:                 'Done',
};

export default defineComponent({
  name: "GeneCallingPage",
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
    Info,
    Dna,
    Download,
    Eye,
    Trash2,
    VueSelect,
    VueSlider,
    Button,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    GeneCallingHelpCollapsible,
  },
  setup() {
    const store = useStore();
    const metag: Ref<boolean> = ref(false);
    const closed_ends: Ref<boolean> = ref(false);
    const mask: Ref<boolean> = ref(false);
    const non_sd: Ref<boolean> = ref(false);
    const min_gene_fraction: Ref<number> = ref(0.10);
    const min_gene_group_fraction: Ref<number> = ref(0.10);
    const tt: Ref<number> = ref(0);
    const numWorkers: Ref<number> = ref(4);
    const uploadedFileNames: Ref<string[]> = ref([]);

    const selectedGenome = ref<string | null>(null);
    const genomeViewerContainer = ref<HTMLDivElement | null>(null);
    let reactRoot: ReturnType<typeof ReactDOM.createRoot> | null = null;
    let activeBlobUrls: string[] = [];
    let blobRangeCache = new Map<string, ArrayBuffer>();
    let fetchPatched = false;

    function revokeBlobUrls() {
        activeBlobUrls.forEach(u => {
            URL.revokeObjectURL(u);
            blobRangeCache.delete(u);
        });
        activeBlobUrls = [];
    }

    function ensureFetchPatchedForBlobRange() {
        if (fetchPatched) return;
        fetchPatched = true;
        const orig = globalThis.fetch.bind(globalThis);
        globalThis.fetch = async (input: RequestInfo | URL, init?: RequestInit) => {
            const url = typeof input === 'string' ? input
                : input instanceof URL ? input.href
                : (input as Request).url;
            if (url.startsWith('blob:')) {
                const hdrs = (init?.headers ?? {}) as Record<string, string>;
                const rangeHdr = hdrs['range'] || hdrs['Range'];
                if (rangeHdr) {
                    const buf = blobRangeCache.get(url);
                    if (buf) {
                        const m = /^bytes=(\d+)-(\d+)?$/.exec(rangeHdr);
                        if (m) {
                            const lo = parseInt(m[1], 10);
                            const hi = m[2] !== undefined ? parseInt(m[2], 10) + 1 : buf.byteLength;
                            const slice = buf.slice(lo, Math.min(hi, buf.byteLength));
                            return new Response(slice, {
                                status: 206,
                                headers: {
                                    'Content-Range': `bytes ${lo}-${Math.min(hi, buf.byteLength) - 1}/${buf.byteLength}`,
                                    'Content-Length': String(slice.byteLength),
                                },
                            });
                        }
                    }
                }
            }
            return orig(input, init);
        };
    }

    function mountGenomeViewer() {
        if (!genomeViewerContainer.value || !selectedGenome.value) return;
        const result = orphosResults.value[selectedGenome.value];
        if (!result) return;

        revokeBlobUrls();

        const fastaUrl = URL.createObjectURL(new Blob([result.fastaBgz]));
        const faiUrl   = URL.createObjectURL(new Blob([result.fastaFai]));
        const gziUrl   = URL.createObjectURL(new Blob([result.fastaGzi]));
        const gffUrl   = URL.createObjectURL(new Blob([result.gffBgz]));
        const csiUrl   = URL.createObjectURL(new Blob([result.gffCsi]));
        activeBlobUrls = [fastaUrl, faiUrl, gziUrl, gffUrl, csiUrl];

        // Pre-populate range-request cache so BgzipFastaAdapter range reads work on blob URLs
        blobRangeCache.set(fastaUrl, result.fastaBgz.buffer as ArrayBuffer);
        blobRangeCache.set(faiUrl,   result.fastaFai.buffer as ArrayBuffer);
        blobRangeCache.set(gziUrl,   result.fastaGzi.buffer as ArrayBuffer);
        blobRangeCache.set(gffUrl,   result.gffBgz.buffer as ArrayBuffer);
        blobRangeCache.set(csiUrl,   result.gffCsi.buffer as ArrayBuffer);
        ensureFetchPatchedForBlobRange();

        const el = React.createElement(GeneCallingGenomeViewer as React.ElementType, {
            fileName: result.fileName,
            fastaUrl,
            faiUrl,
            gziUrl,
            gffUrl,
            csiUrl,
            heightPx: 600,
        });
        if (!reactRoot) {
            reactRoot = ReactDOM.createRoot(genomeViewerContainer.value);
        }
        reactRoot.render(el);
    }

    function selectGenome(fileName: string) {
        selectedGenome.value = fileName;
    }

    onBeforeUnmount(() => {
        if (reactRoot) { reactRoot.unmount(); reactRoot = null; }
        revokeBlobUrls();
    });

    watch(selectedGenome, () => {
        nextTick(() => mountGenomeViewer());
    });

    const { callGenes, resetAllResults_orphos, initCallerWorkers } = useActions(["callGenes", "resetAllResults_orphos", "initCallerWorkers"]);

    function onNumWorkersChange(value: number): void {
      initCallerWorkers(value);
    }

    function resetAll(): void {
      uploadedFileNames.value = [];
      selectedGenome.value = null;
      resetAllResults_orphos();
    }

    function onDropGenome(acceptFiles: File[]): void {
      const newNames = acceptFiles.map(f => f.name);
      for (const name of newNames) {
        if (!uploadedFileNames.value.includes(name)) {
          uploadedFileNames.value.push(name);
        }
      }
      callGenes({
        acceptFiles,
        metag: metag.value,
        closed_ends: closed_ends.value,
        mask: mask.value,
        tt: tt.value,
        non_sd: non_sd.value,
        min_gene_fraction: min_gene_fraction.value,
        min_gene_group_fraction: min_gene_group_fraction.value,
      });
    }

    const orphosResults = computed<Record<string, GeneCallResult>>(() => store.getters.orphosResults);
    const resultCount = computed<number>(() => Object.keys(orphosResults.value).length);

    const geneCallingProgressTotal = computed<number>(() => store.getters.geneCallingProgressTotal);
    const geneCallingStep = computed<string>(() => store.getters.geneCallingStep);
    const geneCallingStepLabel = computed<string>(() => STEP_LABELS[geneCallingStep.value] ?? geneCallingStep.value);

    function downloadGff(result: GeneCallResult): void {
      const blob = new Blob([result.outputFile], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = result.fileName.replace(/\.[^.]+$/, '') + '.gff';
      a.click();
      URL.revokeObjectURL(url);
    }


    function downloadAmrTsv(result: GeneCallResult): void {
      if (!result.amrAnnotationTsv) return;
      const blob = new Blob([result.amrAnnotationTsv], { type: 'text/tab-separated-values' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = result.fileName.replace(/\.[^.]+$/, '') + '.amr.tsv';
      a.click();
      URL.revokeObjectURL(url);
    }

    function downloadZip(): void {
      const files: Record<string, string> = {};
      for (const result of Object.values(orphosResults.value)) {
        const stem = result.fileName.replace(/\.[^.]+$/, '');
        files[stem + '.gff'] = result.outputFile;
        if (result.amrAnnotationTsv) files[stem + '.amr.tsv'] = result.amrAnnotationTsv;
      }
      const bytes = buildZip(files);
      const blob = new Blob([bytes.buffer as ArrayBuffer], { type: 'application/zip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'genecalls.zip';
      a.click();
      URL.revokeObjectURL(url);
    }

    function downloadTarGz(): void {
      const files: Record<string, string> = {};
      for (const result of Object.values(orphosResults.value)) {
        const stem = result.fileName.replace(/\.[^.]+$/, '');
        files[stem + '.gff'] = result.outputFile;
        if (result.amrAnnotationTsv) files[stem + '.amr.tsv'] = result.amrAnnotationTsv;
      }
      const bytes = buildTarGz(files);
      const blob = new Blob([bytes.buffer as ArrayBuffer], { type: 'application/gzip' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'genecalls.tar.gz';
      a.click();
      URL.revokeObjectURL(url);
    }

    const {
      getRootProps: getRootPropsGenome,
      getInputProps: getInputPropsGenome,
      isDragActive: isDragActiveGenome,
      ...restGenome
    } = useDropzone({
      onDrop: onDropGenome,
      accept: fastaExtensionsWithDotAndCompressList,
      multiple: true
    });

    return {
      metag,
      closed_ends,
      mask,
      non_sd,
      tt,
      min_gene_fraction,
      min_gene_group_fraction,
      numWorkers,
      onNumWorkersChange,
      uploadedFileNames,
      resetAll,
      getRootPropsGenome,
      getInputPropsGenome,
      isDragActiveGenome,
      onDropGenome,
      orphosResults,
      resultCount,
      downloadGff,
      downloadAmrTsv,
      downloadZip,
      downloadTarGz,
      selectedGenome,
      genomeViewerContainer,
      selectGenome,
      geneCallingProgressTotal,
      geneCallingStepLabel,
      store,
      ...restGenome,
    };
  },
  computed: {
    orphosError(): string | null {
      return this.store.getters.orphosError;
    },
    genesCalled(): boolean {
      return this.store.getters.genesCalled;
    },
    callingGenes(): boolean {
      return this.store.getters.callingGenes;
    },
    callingGenesFilesSet(): Set<string> {
      return this.store.getters.callingGenesFiles;
    },
  },
  methods: {
    isFileDone(fileName: string): boolean {
      return fileName in (this.orphosResults || {});
    },
    isFileInFlight(fileName: string): boolean {
      return this.callingGenesFilesSet.has(fileName);
    },
  },
});
</script>

<style scoped>
</style>
