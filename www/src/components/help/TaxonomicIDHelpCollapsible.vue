<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Taxonomic ID</span>
      </div>
      <ChevronDown
        class="w-4 h-4 text-gray-500 transition-transform duration-200"
        :class="{ 'rotate-180': isOpen }"
      />
    </CollapsibleTrigger>

    <CollapsibleContent class="px-3 pb-3">
      <Tabs default-value="overview" class="mt-2">
        <TabsList class="grid w-full grid-cols-2 h-8">
          <TabsTrigger value="overview" class="text-xs">Overview</TabsTrigger>
          <TabsTrigger value="parameters" class="text-xs">Parameters</TabsTrigger>
        </TabsList>

        <div class="mt-3 text-sm text-gray-600 max-h-64 overflow-y-auto">
          <TabsContent value="overview" class="space-y-4">
            <p>
              Here you can identify the taxonomic classification at the level of species of bacterial genomic sequences (either
              reads in FASTQ format, or any FASTA sequence). For that, you just need to "upload" (load into memory) your
              sequence files, and you will obtain the three most closest results, with its probability, that the methodology
              obtains for your sequences. You can check the documentation for more information on how this works.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example record and files</p>
              <p class="text-sm">
                <a href="https://www.ebi.ac.uk/ena/browser/view/GCA_004138665.1" target="_blank" class="text-blue-600 hover:underline"><i>Klebsiella pneumoniae</i> assembly: GCA_004138665.1</a>
              </p>
              <ul class="list-disc list-inside mt-1">
                <li>
                  <a href="https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_004138665.1?download=true&gzip=true" class="text-blue-600 hover:underline">GCA_004138665.1.fasta.gz</a>
                </li>
              </ul>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">Minimum Illumina read quality</dt>
                <dd class="ml-4">The quality of the input fastq files, it is an integer from zero to 33. By default, 20.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Minimum counts for k-mer filtering</dt>
                <dd class="ml-4">Another integer for a filtering step, by default 5.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Proportion of reads</dt>
                <dd class="ml-4">When using reads, this allows to subsample them (making thus the analysis faster). A real number between 0 and 1.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Workers</dt>
                <dd class="ml-4">When uploading multiple files, they are processed in different threads. This allows you to select how many you want to create. More will speed up processing of multiple files, but also consume more memory.</dd>
              </div>
            </dl>
          </TabsContent>
        </div>
      </Tabs>
    </CollapsibleContent>
  </Collapsible>
</template>

<script lang="ts">
import { defineComponent, ref } from "vue";
import { HelpCircle, ChevronDown } from "@lucide/vue";
import {
  Collapsible,
  CollapsibleContent,
  CollapsibleTrigger,
} from "@/components/ui/collapsible";
import {
  Tabs,
  TabsContent,
  TabsList,
  TabsTrigger,
} from "@/components/ui/tabs";

export default defineComponent({
  name: "TaxonomicIDHelpCollapsible",
  components: {
    HelpCircle,
    ChevronDown,
    Collapsible,
    CollapsibleContent,
    CollapsibleTrigger,
    Tabs,
    TabsContent,
    TabsList,
    TabsTrigger,
  },
  setup() {
    const isOpen = ref(false);
    return { isOpen };
  },
});
</script>
