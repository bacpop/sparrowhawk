<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium flex items-center gap-2 mb-4">
        <Codesandbox class="w-6 h-6" />
        {{ tabName }}
      </h1>

      <AssemblyHelpCollapsible />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Controls the size of the subsequences (k-mers) used to process the reads.</p>
                </TooltipContent>
              </Tooltip>
              k
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="k"
                         :lazy="true"
                         :min="21"
                         :max="89"
                         :interval="2"
                         :disabled="isProcessingAny"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ k }}
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
                         :disabled="isProcessingAny"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_qual }}
              </span>
            </div>
          </div>

          <div :class="do_fit ? 'opacity-50' : ''">
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
                         :min="infimumMinCount"
                         :max="30"
                         :interval="1"
                         :disabled="isProcessingAny || do_fit"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_count }}
              </span>
            </div>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="do_fit" type="checkbox" v-model="do_fit" :disabled="isProcessingAny"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Allows automatic setting of the minimum k-mer count. Memory usage could increase.</p>
              </TooltipContent>
            </Tooltip>
            <label for="do_fit">
              Automatically set min counts
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="do_bloom" type="checkbox" v-model="do_bloom" :disabled="isProcessingAny"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Alternative filtering with significantly lower memory usage. Recommended for non-small reads. Chunking will be disabled.</p>
              </TooltipContent>
            </Tooltip>
            <label for="do_bloom">
              Use Bloom filter
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="no_deadend" type="checkbox" v-model="no_deadend" :disabled="isProcessingAny"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Maintains short paths in the graph that would usually be pruned.</p>
              </TooltipContent>
            </Tooltip>
            <label for="no_deadend">
              Keep dead-ends
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="no_bubble" type="checkbox" v-model="no_bubble" :disabled="isProcessingAny"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Disables simple bubble collapse in the assembly graph.</p>
              </TooltipContent>
            </Tooltip>
            <label for="no_bubble">
              Do not collapse bubbles
            </label>
          </div>

          <div :class="do_bloom ? 'opacity-50' : ''">
            <label for="csize" class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Size of the preprocessing batches. Lowering this value will usually decrease memory usage somewhat.</p>
                </TooltipContent>
              </Tooltip>
              Chunk size (zero for no chunking)
            </label>
            <input id="csize"
                   type="number"
                   class="w-full border border-gray-300 rounded-md text-sm p-2 active:outline-none focus:outline-none"
                   v-model.number.trim="csize"
                   :disabled="do_bloom || isProcessingAny">
          </div>

          <Button @click="doPreProcess()" v-if="readsProcessed" class="max-w-fit mt-2"
                  variant="outline" size="sm">
            Re-process
          </Button>
        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 flex-1 overflow-hidden pt-12">
      <div v-if="errorInProcessing" class="mx-6 mr-0 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="assemblyErrorType === 'file_count'">
          Assembly requires one or two FASTQ files. Please upload between one and two files and try again.
        </template>
        <template v-else>
          Error found while processing! It is most surely a memory issue: try increasing the chunking, or using a Bloom filter.
        </template>
      </div>

      <!-- Assembly tab -->
      <div v-if="tabName=='Assembly'">

        <!-- Dropzone for file upload - always visible -->
        <div v-bind='getRootPropsReads()'
             :class="[
               'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               isProcessingAny ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer hover:border-gray-400'
             ]">

          <input v-bind='getInputPropsReads()' :disabled="isProcessingAny"/>

          <FileUp/>

          <p v-if='isDragActiveReads'>
            Drop files here ...
          </p>
          <p v-else>
            Drop or click to upload your <b>FASTQ read files</b> (one or two)
          </p>

        </div>

        <!-- Processing status with detailed state -->
        <div v-if="isPreprocessingActive || isAssemblingActive" class="mx-6 mr-0 mt-4 p-3 bg-blue-50 border border-blue-200 rounded-md">
          <div class="flex items-center gap-2 text-sm">
            <Loader2 v-if="isPreprocessingActive" class="w-4 h-4 text-blue-500 animate-spin"/>
            <Loader2 v-else-if="isAssemblingActive" class="w-4 h-4 text-orange-500 animate-spin"/>
            <span class="font-semibold text-gray-800">
              {{ assemblyPhaseTitle }}
            </span>
          </div>
          <div v-if="assemblyStepDetails.step" class="mt-2 ml-6">
            <div class="text-xs font-medium text-gray-700">
              {{ assemblyStepDetails.step }}
            </div>
            <div v-if="assemblyStepDetails.detail" class="text-xs text-gray-600 mt-0.5">
              {{ assemblyStepDetails.detail }}
            </div>
          </div>
        </div>

        <!-- File list with status -->
        <div class="flex flex-row gap-2 w-full mt-4">
          <div v-if="uploadedFiles.length > 0" class="mx-6 w-1/2 flex-grow">
            <div v-for="file in uploadedFiles" :key="file.name"
                 class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">

              <Loader2 v-if="isPreprocessingActive" class="w-4 h-4 text-blue-500 animate-spin"/>
              <Loader2 v-else-if="isAssemblingActive" class="w-4 h-4 text-orange-500 animate-spin"/>
              <Check v-else-if="readsProcessed || readsPreprocessed" class="w-4 h-4 text-green-500"/>

              <span class="flex-grow text-sm font-mono truncate">
              {{ file.name }}
            </span>
            </div>
          </div>
          <div v-if="readsProcessed" class="w-1/2 flex-grow">
            <DownloadButton/>
          </div>
        </div>

        <!-- Reset button when processing is done or in progress -->
        <div v-if="uploadedFiles.length > 0" class="mx-6 mt-4">
          <Button variant="outline" size="sm" @click="resetAll">
            <Trash2 class="mr-1 h-3 w-3" />
            Clear results
          </Button>
        </div>

        <slot/>
      </div>
    </div>
  </div>
</template>


<script lang="ts">
import {defineComponent, ref, Ref} from "vue";
import {useDropzone} from "vue3-dropzone";
import {useActions, useState} from "vuex-composition-helpers";
import {useStore} from "vuex";
import VueSlider from 'vue-3-slider-component';
import "@fontsource/ibm-plex-mono";
import {Check, FileUp, Loader2, Info, Codesandbox, Trash2} from "lucide-vue-next";
import {Tooltip, TooltipContent, TooltipProvider, TooltipTrigger} from "@/components/ui/tooltip";
import DownloadButton from "@/components/DownloadButton.vue";
import AssemblyHelpCollapsible from "@/components/help/AssemblyHelpCollapsible.vue";
import {Button} from "@/components/ui/button";
import { fastqExtensionsWithDotAndCompressList } from "@/utils";

export default defineComponent({
  name: "AssemblyPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    DownloadButton,
    AssemblyHelpCollapsible,
    Button,
    VueSlider,
    FileUp,
    Loader2,
    Check,
    Info,
    Codesandbox,
    Trash2,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger
  },
  setup() {
    const store = useStore();
    const k: Ref<number> = ref(31);
    const min_count: Ref<number> = ref(5);
    const min_qual: Ref<number> = ref(20);
    const csize: Ref<number> = ref(150000);
    const do_bloom: Ref<boolean> = ref(false);
    const do_fit: Ref<boolean> = ref(true);
    const assemblying: Ref<boolean> = ref(false);
    const no_bubble: Ref<boolean> = ref(false);
    const no_deadend: Ref<boolean> = ref(false);
    const uploadedFiles: Ref<File[]> = ref([]);

    const {
      processReads,
      doTheAssembly,
      resetAllResults,
      removeErrors
    } = useActions(["processReads", "doTheAssembly", "resetAllResults", "removeErrors"]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const {allResults} = useState(["allResults"]) as any;

    function onDropReads(acceptFiles: File[]): void {
      uploadedFiles.value = acceptFiles;
      doPreProcess();
    }

    function doPreProcess(): void {
      processReads({
        acceptFiles: uploadedFiles.value,
        k: k.value,
        min_count: min_count.value,
        min_qual: min_qual.value,
        csize: csize.value,
        do_bloom: do_bloom.value,
        do_fit: do_fit.value,
        no_bubble_collapse: no_bubble.value,
        no_dead_end_removal: no_deadend.value
      });
    }

    function doAss(): void {
      assemblying.value = true;
      doTheAssembly();
    }

    function resetAll(): void {
      assemblying.value = false;
      uploadedFiles.value = [];
      resetAllResults();
      removeErrors();
    }

    const {
      getRootProps: getRootPropsReads,
      getInputProps: getInputPropsReads,
      isDragActive: isDragActiveReads,
      ...restReads
    } = useDropzone({
      onDrop: onDropReads,
      accept: fastqExtensionsWithDotAndCompressList,
      multiple: true
    });

    return {
      store,
      k,
      min_count,
      min_qual,
      csize,
      do_bloom,
      do_fit,
      no_bubble,
      no_deadend,
      assemblying,
      resetAll,
      doPreProcess,
      doAss,
      getRootPropsReads,
      getInputPropsReads,
      isDragActiveReads,
      onDropReads,
      allResults,
      uploadedFiles,
      ...restReads,
    };
  },

  computed: {
    readsProcessed(): boolean {
      return this.store.getters.queryAssembled;
    },
    readsProcessing(): boolean {
      return this.store.getters.readsProcessing;
    },
    readsPreprocessed(): boolean {
      return this.store.getters.readsPreprocessed;
    },
    readsPreprocessing(): boolean {
      return this.store.getters.readsPreprocessing;
    },
    errorInProcessing(): boolean {
      return this.store.getters.getErrors;
    },
    assemblyErrorType(): string {
      return this.store.getters.getErrorType;
    },
    readsName(): string {
      return this.store.getters.readsName;
    },
    isPreprocessingActive(): boolean {
      return this.store.getters.isPreprocessing;
    },
    isAssemblingActive(): boolean {
      return this.store.getters.isAssembling;
    },
    isProcessingAny(): boolean {
      return this.store.getters.isPreprocessing || this.store.getters.isAssembling;
    },
    infimumMinCount() {
      return this.do_bloom ? 3 : 0
    },
    currentAssemblyState(): string {
      return this.store.getters.assemblyState;
    },
    assemblyPhaseTitle(): string {
      if (this.isPreprocessingActive) return 'Pre-processing';
      if (this.isAssemblingActive) return 'Assembling';
      return '';
    },
    assemblyStepDetails(): { step: string; detail: string } {
      const state = this.currentAssemblyState;
      if (!state) return { step: '', detail: '' };

      // Parse preprocessing states
      if (state.startsWith('preprocess:')) {
        // Initialization
        if (state === 'preprocess:start') {
          return { step: 'Initialising', detail: 'Setting up preprocessing pipeline' };
        }

        // Bloom filter preprocessing
        if (state.includes('bloom')) {
          if (state === 'preprocess:bloom:start') {
            return { step: 'Bloom Filter Setup', detail: 'Initialising Bloom filter' };
          }
          if (state === 'preprocess:bloom:loop:start') {
            return { step: 'Reading FASTQ Files', detail: 'Processing paired-end reads with Bloom filter' };
          }
          
          // Progress messages
          if (state.includes('preprocess:bloom:loop:')) {
            const parts = state.split(':');
            if (parts.length >= 4) {
              const numReads = isNaN(Number(parts[3])) ? 0 : Number(parts[3]);
              const percentage = parts[4];
              
              if (percentage) {
                return { 
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads (${percentage}% complete)`
                };
              } else {
                return { 
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads`
                };
              }
            }
          }
          
          if (state === 'preprocess:bloom:loop:end') {
            return { step: 'Reading FASTQ Files', detail: 'All reads processed successfully' };
          }
          if (state === 'preprocess:bloom:fitting') {
            return { step: 'Optimising Parameters', detail: 'Fitting k-mer count threshold' };
          }
          if (state === 'preprocess:bloom:filtering') {
            return { step: 'Filtering K-mers', detail: 'Filtering k-mers from dataset' };
          }
        }
        
        // Bulk (no chunk, no bloom) preprocessing
        else if (state.includes('bulk')) {
          if (state === 'preprocess:bulk:start') {
            return { step: 'Bulk Processing Setup', detail: 'Initialising' };
          }
          if (state === 'preprocess:bulk:loop:start') {
            return { step: 'Reading FASTQ Files', detail: 'Processing paired-end reads' };
          }
          
          // Progress messages
          if (state.includes('preprocess:bulk:loop:')) {
            const parts = state.split(':');
            if (parts.length >= 4) {
              const numReads = (isNaN(Number(parts[3])) ? 0 : Number(parts[3]));
              const percentage = parts[4];
              
              if (percentage) {
                return { 
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads (${percentage}% complete)`
                };
              } else {
                return { 
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads`
                };
              }
            }
          }
          
          if (state === 'preprocess:bulk:loop:end') {
            return { step: 'Reading FASTQ Files', detail: 'All reads processed successfully' };
          }
          if (state === 'preprocess:bulk:sorting') {
            return { step: 'Sorting K-mers', detail: 'Organising k-mers for efficient filtering' };
          }
          if (state === 'preprocess:bulk:fitting') {
            return { step: 'Optimising Parameters', detail: 'Automatically fitting k-mer count threshold' };
          }
          if (state === 'preprocess:bulk:filtering') {
            return { step: 'Filtering K-mers', detail: 'Removing low-quality k-mers from dataset' };
          }
        }
        
        // Chunked preprocessing
        else if (state.includes('chunked')) {
          if (state === 'preprocess:chunked:start') {
            return { step: 'Chunked Processing Setup', detail: 'Initialising memory-optimised chunked mode' };
          }
          if (state === 'preprocess:chunked:fitting') {
            return { step: 'Optimising Parameters', detail: 'Automatically fitting k-mer count threshold across chunks' };
          }
          if (state === 'preprocess:chunked:filtering') {
            return { step: 'Filtering K-mers', detail: 'Removing low-quality k-mers from chunked dataset' };
          }

          // Progress messages
          if (state.includes('preprocess:chunked:loop:')) {
            const parts = state.split(':');
            if (parts.length >= 4) {
              const numReads = isNaN(Number(parts[3])) ? 0 : Number(parts[3]);
              const percentage = parts[4];

              if (percentage) {
                return {
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads (${percentage}% complete)`
                };
              } else {
                return {
                  step: 'Reading FASTQ Files',
                  detail: `Processed ${numReads.toLocaleString()} reads`
                };
              }
            }
          }
        }
        
        // Common preprocessing states
        if (state === 'preprocess:saving') {
          return { step: 'Saving Data', detail: 'Writing preprocessed k-mers to memory' };
        }
        if (state === 'preprocess:end') {
          return { step: 'Complete', detail: 'Preprocessing finished successfully' };
        }
      }
      
      // Parse assembly states
      if (state.startsWith('assembly:')) {
        if (state === 'assembly:start') {
          return { step: 'Starting Assembly', detail: 'Initialising assembly pipeline' };
        }
        if (state === 'assembly:create_graph') {
          return { step: 'Building Graph', detail: 'Creating de Bruijn graph' };
        }
        if (state === 'assembly:correct_graph') {
          return { step: 'Correcting Graph', detail: 'Processing' };
        }
        if (state === 'assembly:collapse_graph') {
          return { step: 'Collapsing Graph', detail: 'Processing' };
        }
        if (state === 'assembly:saving') {
          return { step: 'Saving Results', detail: 'Saving to FASTA format' };
        }
        if (state === 'assembly:end') {
          return { step: 'Complete', detail: 'Assembly finished successfully' };
        }
      }
      
      // Initialization
      if (state === 'initialised') {
        return { step: 'Initialised', detail: 'System ready' };
      }
      
      // Fallback to raw state if no match
      return { step: state, detail: '' };
    }
  },

  watch: {
    readsPreprocessed(newVal: boolean) {
      if (newVal) {
        this.doAss();
      }
    },
    do_bloom(newVal: boolean) {
      if (newVal && this.min_count < 3) {
        this.min_count = 3;
      }
    }
  },

  methods: {
    clear(): void {
      this.resetAll();
    },

    doAssembly(): void {
      this.doAss();
    },
  },
});

</script>

<style scoped>
</style>
