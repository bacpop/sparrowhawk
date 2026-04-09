<template>
  <div class="flex flex-row">
    <div class="w-[350px] shrink-0">
      <h1 class="text-2xl font-medium mb-4 flex items-center gap-2">
        <TextAlignCenter v-if="tabName === 'Mapping'" class="w-6 h-6" />
        <TreePine v-else-if="tabName === 'Alignment'" class="w-6 h-6" />
        <div v-if="tabName === 'Mapping'">
            Mapping (within species)
        </div>
        <div v-if="tabName === 'Alignment'">
            Alignment (within species)
        </div>
      </h1>

      <MappingHelpCollapsible v-if="tabName === 'Mapping'" />
      <AlignmentHelpCollapsible v-else-if="tabName === 'Alignment'" />

      <TooltipProvider>
        <div class="flex flex-col gap-4">
          <div>
            <p class="flex items-center gap-1">
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Controls the size of the subsequences (k-mers) used to process the reads. Maximum value is 63 for this mode. You can only set it before uploading any file.</p>
                </TooltipContent>
              </Tooltip>
              k
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSlider class="flex-grow"
                         v-model="k"
                         :lazy="true"
                         :min="5"
                         :max="63"
                         :interval="2"
                         :disabled="refProcessed || hasAlignmentResults || isProcessingAny"
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
          
          <!-- TODO: implement properly minimum count for allele frequency, like in Ska you can do -->
          
          <!-- <div>
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
                         :min="0"
                         :max="30"
                         :interval="1"
                         :disabled="isProcessingAny"
              />
              <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                {{ min_count }}
              </span>
            </div>
          </div> -->

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
                         :disabled="isProcessingAny"
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
                  <p class="max-w-xs">Either no filtering, filter only the middle base of split k-mers, or any nucleotide in a k-mer.</p>
                </TooltipContent>
              </Tooltip>
              Quality filter type
            </p>
            <div class="flex flex-row items-center w-full gap-2">
              <VueSelect class="flex-grow"
                         v-model="qual_filter"
                         :options="[
                          { label: 'No filter', value: 0},
                          { label: 'Middle base', value: 1},
                          { label: 'All bases', value: 2},
                         ]"
                         :isDisabled="isProcessingAny"
                         :isClearable=false
              />
            </div>
          </div>

          <div class="flex flex-row items-center w-full gap-2">
            <input id="rc" type="checkbox" v-model="rc" :disabled="refProcessed || hasAlignmentResults"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Recommended to activate when using reads. You can only set it before uploading any file.</p>
              </TooltipContent>
            </Tooltip>
            <label for="rc">
              Use canonical k-mers
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2" v-if="tabName=='Mapping'">
            <input id="ambig_mask" type="checkbox" v-model="ambig_mask" :disabled="refProcessed || hasAlignmentResults"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Replace amiguous bases in the alignment with "N".</p>
              </TooltipContent>
            </Tooltip>
            <label for="ambig_mask">
              Mask ambiguous bases
            </label>
          </div>

          <div class="flex flex-row items-center w-full gap-2" v-if="tabName=='Mapping'">
            <input id="repeat_mask" type="checkbox" v-model="repeat_mask" :disabled="refProcessed || hasAlignmentResults"/>
            <Tooltip>
              <TooltipTrigger as-child>
                <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
              </TooltipTrigger>
              <TooltipContent>
                <p class="max-w-xs">Mask any repeats in the alignment with "N".</p>
              </TooltipContent>
            </Tooltip>
            <label for="repeat_mask">
              Mask repeats
            </label>
          </div>

          <!-- Transmission clustering options (Alignment tab only) -->
          <template v-if="tabName === 'Alignment'">
            <hr class="border-gray-200" />

            <div class="flex flex-row items-center w-full gap-2">
              <input id="enable_clustering" type="checkbox" v-model="enableClustering" :disabled="isProcessingAny || isClustering"/>
              <Tooltip>
                <TooltipTrigger as-child>
                  <Info class="w-3.5 h-3.5 text-gray-400 cursor-help" />
                </TooltipTrigger>
                <TooltipContent>
                  <p class="max-w-xs">Perform single-linkage transmission clustering on the alignment and display an interactive transmission graph.</p>
                </TooltipContent>
              </Tooltip>
              <label for="enable_clustering">
                Transmission clustering
              </label>
            </div>

            <div v-if="enableClustering">
              <p class="flex items-center gap-1">
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
                <VueSlider class="flex-grow"
                           v-model="snp_threshold"
                           :lazy="true"
                           :min="1"
                           :max="1000"
                           :interval="1"
                           :disabled="isClustering"
                />
                <span class="block w-[40px] text-center border border-gray-300 rounded-md text-sm">
                  {{ snp_threshold }}
                </span>
              </div>

              <Button
                v-if="hasAlignmentResults && !isClustering"
                @click="runClustering"
                class="mt-2 w-full"
                variant="outline"
                size="sm"
              >
                <Network class="mr-1 h-3 w-3" />
                Run transmission clustering
              </Button>
              <div v-else-if="isClustering" class="mt-2 flex items-center gap-2 text-sm text-gray-500">
                <Loader2 class="w-4 h-4 animate-spin text-amber-500" />
                Clustering...
              </div>

              <hr class="border-gray-200 mt-3" />

              <!-- Metadata upload -->
              <div class="mt-2">
                <p class="text-sm mb-1">Optional metadata file (CSV)</p>
                <label class="block cursor-pointer">
                  <span class="block p-2 border border-gray-200 rounded-md text-sm text-gray-500 hover:border-gray-400 text-center">
                    {{ metadataFileName || 'Upload metadata CSV' }}
                  </span>
                  <input type="file" accept=".csv" class="hidden" @change="onMetadataFileChange" />
                </label>
                <div v-if="metadataError"
                     class="mt-2 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
                  {{ metadataError }}
                </div>
              </div>

              <!-- Location matching -->
              <div v-if="metadataRows.length > 0" class="mt-2 flex items-center gap-2">
                <input id="enable_location" type="checkbox" v-model="enableLocationMatching" />
                <label for="enable_location" class="text-sm">Enable location matching</label>
              </div>
            </div>
          </template>

        </div>
      </TooltipProvider>
    </div>

    <div class="min-w-0 flex-1 overflow-hidden pt-12">

      <div v-if="skaError"
           class="mx-6 mb-4 p-3 bg-red-50 border border-red-300 rounded-md text-sm text-red-800">
        <template v-if="skaError === 'memory'">
          Error during processing — most likely a memory issue. Try with fewer or smaller files.
        </template>
        <template v-else>
          An unexpected error occurred. Please reset and try again.
        </template>
      </div>

      <!-- Mapping tab -->
      <div v-if="tabName=='Mapping'">

        <!-- Single dropzone - always visible -->
        <div v-bind='getRootPropsMapping()'
             :class="[
               'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               isProcessingAny ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind='getInputPropsMapping()' :disabled="isProcessingAny"/>
          <FileUp/>
          <p v-if='isDragActiveMapping'>
            Drop files here ...
          </p>
          <p v-else>
            Drop or click to upload your <b>{{ dropzonePrompt }}</b>
          </p>
        </div>

        <p v-if="isIndexingRef" class="mx-6 mr-0 mt-4 text-sm text-gray-500">
          Indexing reference...
        </p>
        <p v-else-if="isMapping" class="mx-6 mr-0 mt-4 text-sm text-gray-500">
          Mapping...
        </p>

        <!-- File list with status and download button -->
        <div class="flex flex-row gap-2 w-full mt-4">
          <div v-if="uploadedFiles.length > 0" class="mx-6 w-1/2 flex-grow">
            <div v-for="file in uploadedFiles" :key="file.name"
                 class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">

              <Loader2 v-if="getFileStatus(file) === 'indexing'" class="w-4 h-4 text-blue-500 animate-spin"/>
              <Loader2 v-else-if="getFileStatus(file) === 'mapping'" class="w-4 h-4 text-orange-500 animate-spin"/>
              <Check v-else-if="getFileStatus(file) === 'done'" class="w-4 h-4 text-green-500"/>

              <span class="flex-grow text-sm font-mono truncate">
                {{ file.name }}
              </span>
              <span v-if="file.type === 'reference'" class="text-xs text-gray-400">(ref)</span>
            </div>
          </div>
          <div v-if="hasMappingResults" class="w-1/2 flex-grow">
            <DownloadButtonSka />
          </div>
        </div>
        

        <!-- Reset button -->
        <Button v-if="uploadedFiles.length > 0" @click="resetAll" class="mx-6 mr-0 mt-4" variant="outline" size="sm">
          <Trash2 class="mr-1 h-3 w-3" />
          Clear results
        </Button>

        <div v-if="hasMappingResults" class="mx-6 mt-6">
          <div class="w-full">
            <MSAViewer
              :data="mappingMSAData"
              filename="Mapping Results"
              :show-header="false"
              class="h-[400px] w-full"
            />
          </div>
        </div>

        <slot name="mapping"/>
      </div>

      <!-- Alignment tab -->
      <div v-else-if="tabName=='Alignment'">
        <!-- Dropbox - always visible when not aligning -->
        <div v-if="!isAligning"
             v-bind='getRootPropsQueryAlign()'
             :class="[
               'p-6 mx-6 mr-0 bg-white border border-gray-200 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600',
               'cursor-pointer hover:border-gray-400'
             ]">
          <input v-bind='getInputPropsQueryAlign()'/>
          <FileUp/>
          <p v-if='isDragActiveQueryAlign'>
            Drop files here ...
          </p>
          <p v-else>
            Drop or click to upload <b>files to be aligned</b>
          </p>
        </div>

        <div v-else class="p-6 mx-6 bg-amber-50 border border-amber-400 rounded-md flex flex-col justify-center items-center gap-2 text-gray-600">
          <Loader2 class="w-6 h-6 text-amber-500 animate-spin"/>
          <p class="text-sm text-gray-500">Aligning...</p>
        </div>

        <!-- File list with status -->
        <div v-if="uploadedAlignmentFiles.length > 0" class="mx-6 mr-0 mt-4 max-h-48 overflow-y-auto">
          <div v-for="fileName in uploadedAlignmentFiles" :key="fileName"
               class="flex items-center gap-2 py-2 px-3 bg-gray-50 rounded-md mb-2">
            <Loader2 v-if="isAligning" class="w-4 h-4 text-orange-500 animate-spin"/>
            <Check v-else-if="hasAlignmentResults" class="w-4 h-4 text-green-500"/>
            <span class="flex-grow text-sm font-mono truncate">
              {{ fileName }}
            </span>
          </div>
        </div>


        <div v-if="hasAlignmentResults" class="px-8">
          <DownloadButtonSkaAlignment />
        </div>

        <!-- Reset button -->
        <Button v-if="uploadedAlignmentFiles.length > 0 && !isProcessingAny" @click="resetAll" class="mx-6 mt-4" variant="outline" size="sm">
          <Trash2 class="mr-1 h-3 w-3" />
          Clear results
        </Button>

        <!-- Transmission cluster results -->
        <div v-if="hasClusterResults" class="mx-6 mt-6">
          <div class="flex items-center justify-between mb-2">
            <h2 class="text-base font-medium">Transmission Clusters</h2>
            <Button @click="downloadClustersTSV" variant="outline" size="sm">
              <FileDown class="mr-1 h-3 w-3" />
              Download TSV
            </Button>
          </div>

          <div class="max-h-48 overflow-y-auto border border-gray-200 rounded-md text-sm">
            <table class="w-full">
              <thead class="bg-gray-50 sticky top-0">
                <tr>
                  <th class="text-left px-3 py-2 font-medium text-gray-600">Sample</th>
                  <th class="text-left px-3 py-2 font-medium text-gray-600">Cluster</th>
                </tr>
              </thead>
              <tbody>
                <tr v-for="row in clusterTableRows" :key="row.name"
                    class="border-t border-gray-100 hover:bg-gray-50">
                  <td class="px-3 py-2 font-mono truncate max-w-xs">{{ row.name }}</td>
                  <td class="px-3 py-2">{{ row.cluster }}</td>
                </tr>
              </tbody>
            </table>
          </div>

          <div class="mt-4 h-[500px]">
            <TransmissionGraph
              :nodes="transmissionGraphNodes"
              :links="transmissionGraphLinks"
            />
          </div>

          <div v-if="metadataRows.length > 0" class="mt-6">
            <h3 class="text-sm font-medium text-gray-600 mb-2">Cluster timeline</h3>
            <TransmissionClusterTimeline
              :clusterRows="clusterTableRows"
              :metadataRows="metadataRows"
            />
          </div>
        </div>

        <slot name="alignment"/>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, Ref, computed } from "vue";
import { useDropzone } from "vue3-dropzone";
import { useActions, useState } from "vuex-composition-helpers";
import { useStore } from "vuex";
import VueSlider from 'vue-3-slider-component';
import VueSelect from "vue3-select-component";
import "vue3-select-component/styles";
import { Check, FileUp, FileDown, Loader2, Info, Network, TextAlignCenter, TreePine, Trash2 } from "lucide-vue-next";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Button } from "@/components/ui/button";
import MappingHelpCollapsible from "@/components/help/MappingHelpCollapsible.vue";
import AlignmentHelpCollapsible from "@/components/help/AlignmentHelpCollapsible.vue";
import DownloadButtonSka from "@/components/SequenceViewer/DownloadButtonSka.vue";
import DownloadButtonSkaAlignment from "@/components/SequenceViewer/DownloadButtonSkaAlignment.vue";
import { MSAViewer } from "@/components/MSAViewer";
import TransmissionGraph from "@/components/TransmissionGraph.vue";
import TransmissionClusterTimeline from "@/components/TransmissionClusterTimeline.vue";
import { fastxExtensionsWithDotAndCompressList } from "@/utils";

interface UploadedFile {
  name: string;
  type: 'reference' | 'query';
}

export default defineComponent({
  name: "MappingAlignmentPage",
  props: {
    tabName: {
      type: String,
      required: true
    }
  },
  components: {
    VueSlider,
    VueSelect,
    FileUp,
    FileDown,
    Loader2,
    Check,
    Info,
    Network,
    TextAlignCenter,
    TreePine,
    Trash2,
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
    Button,
    MappingHelpCollapsible,
    AlignmentHelpCollapsible,
    DownloadButtonSka,
    DownloadButtonSkaAlignment,
    MSAViewer,
    TransmissionGraph,
    TransmissionClusterTimeline,
  },
  setup() {
    const store = useStore();
    const k: Ref<number> = ref(31);
    const min_qual: Ref<number> = ref(20);
    const min_count: Ref<number> = ref(0);
    const proportion_reads: Ref<number> = ref(1);
    const qual_filter: Ref<number> = ref(2);
    const rc: Ref<boolean> = ref(false);
    const ambig_mask: Ref<boolean> = ref(false);
    const repeat_mask: Ref<boolean> = ref(false);
    const enableClustering: Ref<boolean> = ref(false);
    const snp_threshold: Ref<number> = ref(20);
    const uploadedFiles: Ref<UploadedFile[]> = ref([]);
    const uploadedAlignmentFiles: Ref<string[]> = ref([]);
    const metadataRows: Ref<{ id: string; date: string; location: string }[]> = ref([]);
    const metadataError: Ref<string | null> = ref(null);
    const metadataFileName: Ref<string | null> = ref(null);
    const enableLocationMatching: Ref<boolean> = ref(false);


    const {
      processRef,
      processQueryMap,
      processQueryAlign,
      processCluster,
      resetAllResults_ska
    } = useActions(["processRef", "processQueryMap", "processQueryAlign", "processCluster", "resetAllResults_ska"]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const { allResults_ska } = useState(["allResults_ska"]) as any;

    const refProcessed = computed(() => store.getters.refProcessed);
    const isIndexingRef = computed(() => store.getters.isIndexingRef);
    const isMapping = computed(() => store.getters.isMapping);

    function onDropMapping(acceptFiles: File[]): void {
      if (!refProcessed.value && !isIndexingRef.value) {
        // First upload is reference
        uploadedFiles.value = [{ name: acceptFiles[0].name, type: 'reference' }];
        processRef({
          acceptFiles: acceptFiles,
          k: k.value, 
          rc: rc.value,
          ambig_mask: ambig_mask.value,
          repeat_mask: repeat_mask.value,
        });
      } else if (refProcessed.value && !isMapping.value) {
        // Subsequent uploads are query files
        const newFiles = acceptFiles.map(f => ({ name: f.name, type: 'query' as const }));
        uploadedFiles.value = [...uploadedFiles.value, ...newFiles];
        processQueryMap({ 
          acceptFiles: acceptFiles, 
          proportion_reads: proportion_reads.value,
          min_count: min_count.value,
          min_qual: min_qual.value,
          qual_filter: qual_filter.value, });
      }
    }

    function onDropQueryAlign(acceptFiles: File[]): void {
      const newNames = acceptFiles.map(f => f.name);
      // Append new files to the existing list, avoiding duplicates
      const existingNames = new Set(uploadedAlignmentFiles.value);
      for (const name of newNames) {
        if (!existingNames.has(name)) {
          uploadedAlignmentFiles.value.push(name);
        }
      }
      processQueryAlign({ 
        acceptFiles: acceptFiles,
        k: k.value,
        proportion_reads: proportion_reads.value, 
        rc: rc.value,
        min_count: min_count.value,
        min_qual: min_qual.value,
        qual_filter: qual_filter.value
     });
    }

    function parseMetadataCSV(file: File): void {
      metadataError.value = null;
      metadataFileName.value = file.name;
      const reader = new FileReader();
      reader.onload = (e) => {
        const text = e.target?.result as string;
        const lines = text.split('\n').filter(l => l.trim());
        if (lines.length < 2) { metadataError.value = "Metadata file is empty or has no data rows."; return; }
        const headers = lines[0].split(',').map(h => h.trim().toLowerCase());
        const idIdx = headers.indexOf('id');
        const dateIdx = headers.indexOf('date');
        const locIdx = headers.indexOf('location');
        if (idIdx === -1 || dateIdx === -1 || locIdx === -1) {
          metadataError.value = "Metadata CSV must have columns: ID, date, location.";
          return;
        }
        const datePattern = /^\d{4}-\d{2}-\d{2}$/;
        const rows: { id: string; date: string; location: string }[] = [];
        const errors: string[] = [];
        for (let i = 1; i < lines.length; i++) {
          const cols = lines[i].split(',');
          const id = cols[idIdx]?.trim();
          const date = cols[dateIdx]?.trim() ?? '';
          const location = cols[locIdx]?.trim() ?? '';
          if (!id) continue;
          if (date && !datePattern.test(date)) {
            errors.push(`Row ${i + 1}: invalid date "${date}" (expected yyyy-mm-dd)`);
          }
          rows.push({ id, date, location });
        }
        if (errors.length > 0) {
          metadataError.value = errors.join('; ');
          return;
        }
        metadataRows.value = rows;
      };
      reader.readAsText(file);
    }

    function onMetadataFileChange(event: Event): void {
      const input = event.target as HTMLInputElement;
      if (input.files && input.files[0]) {
        parseMetadataCSV(input.files[0]);
      }
    }

    function resetAll(): void {
      uploadedFiles.value = [];
      uploadedAlignmentFiles.value = [];
      enableClustering.value = false;
      snp_threshold.value = 20;
      metadataRows.value = [];
      metadataError.value = null;
      metadataFileName.value = null;
      enableLocationMatching.value = false;
      resetAllResults_ska();
    }

    const {
      getRootProps: getRootPropsMapping,
      getInputProps: getInputPropsMapping,
      isDragActive: isDragActiveMapping,
      ...restMapping
    } = useDropzone({
      onDrop: onDropMapping,
      accept: fastxExtensionsWithDotAndCompressList
    });

    const {
      getRootProps: getRootPropsQueryAlign,
      getInputProps: getInputPropsQueryAlign,
      isDragActive: isDragActiveQueryAlign,
      ...restQueryAlign
    } = useDropzone({
      onDrop: onDropQueryAlign,
      accept: fastxExtensionsWithDotAndCompressList
    });

    return {
      store,
      k,
      min_qual,
      min_count,
      proportion_reads,
      qual_filter,
      rc,
      ambig_mask,
      repeat_mask,
      enableClustering,
      snp_threshold,
      resetAll,
      uploadedFiles,
      uploadedAlignmentFiles,
      processCluster,
      getRootPropsMapping,
      getInputPropsMapping,
      isDragActiveMapping,
      getRootPropsQueryAlign,
      getInputPropsQueryAlign,
      isDragActiveQueryAlign,
      allResults_ska,
      metadataRows,
      metadataError,
      metadataFileName,
      enableLocationMatching,
      onMetadataFileChange,
      ...restMapping,
      ...restQueryAlign
    };
  },
  computed: {
    skaError(): string | null {
      return this.store.getters.skaError;
    },
    refProcessed(): boolean {
      return this.store.getters.refProcessed;
    },
    refName(): string {
      return this.store.getters.refName;
    },
    isIndexingRef(): boolean {
      return this.store.getters.isIndexingRef;
    },
    isMapping(): boolean {
      return this.store.getters.isMapping;
    },
    isAligning(): boolean {
      return this.store.getters.isAligning;
    },
    isProcessingAny(): boolean {
      return this.store.getters.isIndexingRef || this.store.getters.isMapping || this.store.getters.isAligning;
    },
    isClustering(): boolean {
      return this.store.getters.isClustering;
    },
    hasClusterResults(): boolean {
      return this.store.getters.hasClusterResults;
    },
    adjustedClusterResults(): Record<string, number> | null {
      const original = this.allResults_ska.clusterResults;
      if (!original || !this.enableLocationMatching || !this.metadataRows.length) return original;
      const keyToId = new Map<string, number>();
      let nextId = 1;
      const result: Record<string, number> = {};
      for (const [name, cluster] of Object.entries(original)) {
        const loc = this.sampleLocation(name) ?? '';
        const key = `${cluster}:${loc}`;
        if (!keyToId.has(key)) keyToId.set(key, nextId++);
        result[name] = keyToId.get(key)!;
      }
      return result;
    },
    clusterTableRows(): { name: string; cluster: number }[] {
      const results = this.adjustedClusterResults;
      if (!results) return [];
      return Object.entries(results)
        .map(([name, cluster]) => ({ name, cluster: cluster as number }))
        .sort((a, b) => a.cluster - b.cluster || a.name.localeCompare(b.name));
    },
    transmissionGraphNodes() {
      const nodes = this.allResults_ska.transmissionGraph?.nodes ?? [];
      const adjusted = this.adjustedClusterResults;
      if (!adjusted) return nodes;
      return nodes.map((n: { id: string; cluster: number }) => ({
        ...n,
        cluster: adjusted[n.id] ?? n.cluster,
      }));
    },
    transmissionGraphLinks() {
      return this.allResults_ska.transmissionGraph?.links ?? [];
    },
    hasAlignmentResults(): boolean {
      return this.allResults_ska.alignResults[0] && this.allResults_ska.alignResults[0].aligned;
    },
    hasMappingResults(): boolean {
      return Object.keys(this.allResults_ska.mapResults).length > 0;
    },
    // Convert mapping results to MSA format for the viewer
    mappingMSAData(): { id: string; sequence: string }[] {
      const results: { id: string; sequence: string }[] = [];

      // Add reference sequence first
      if (this.allResults_ska.ref && this.allResults_ska.ref.length > 0) {
        const refSeq = this.allResults_ska.ref.join('');
        console.log('[mappingMSAData] ref sequence length:', refSeq.length);
        results.push({
          id: this.refName || 'Reference',
          sequence: refSeq
        });
      }

      // Add mapped sequences
      for (const fileName of Object.keys(this.allResults_ska.mapResults)) {
        const mapping = this.allResults_ska.mapResults[fileName];
        if (mapping?.mapped_sequences && mapping.mapped_sequences.length > 0) {
          const seq = mapping.mapped_sequences.join('');
          console.log('[mappingMSAData] mapped sequence length:', seq.length, 'for', fileName);
          results.push({
            id: fileName,
            sequence: seq
          });
        }
      }

      console.log('[mappingMSAData] total sequences:', results.length);
      return results;
    },
    dropzonePrompt(): string {
      if (!this.refProcessed && !this.isIndexingRef) {
        return 'reference fasta file';
      }
      return 'files to be mapped';
    }
  },

  methods: {
    clear(): void {
      this.resetAll();
    },
    sampleLocation(name: string): string | null {
      const prefix = name.replace(/(?:\.[^.]+)+$/, '');
      const row = (this.metadataRows as { id: string; date: string; location: string }[])
        .find(r => r.id === prefix || r.id === name);
      return row?.location ?? null;
    },
    runClustering(): void {
      this.processCluster({ snp_threshold: this.snp_threshold });
    },
    downloadClustersTSV(): void {
      const rows = this.clusterTableRows;
      const tsv = "Sample\tCluster\n" + rows.map(r => `${r.name}\t${r.cluster}`).join("\n");
      const blob = new Blob([tsv], { type: "text/tab-separated-values" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "transmission_clusters.tsv";
      a.click();
      URL.revokeObjectURL(url);
    },
    getFileStatus(file: UploadedFile): 'indexing' | 'mapping' | 'done' {
      if (file.type === 'reference') {
        if (this.isIndexingRef) return 'indexing';
        if (this.refProcessed) return 'done';
      } else {
        // Query file - check if mapped_sequences has data
        const result = this.allResults_ska.mapResults[file.name];
        if (this.isMapping) return 'mapping';
        if (result?.mapped_sequences?.length > 0 || !this.isMapping) return 'done';
      }
      return 'mapping';
    }
  },
});
</script>

<style scoped>
</style>
