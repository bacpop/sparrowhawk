<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Host depletion</span>
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
              Here you can either select or filter your reads (or in general, sequences) to, for example, remove host reads from your sequencing experiment data.
              To do so, you will need an index, to be created with <a href="https://github.com/bede/deacon" target="_blank" class="text-blue-600 hover:underline">Deacon</a>. You have an example index below you can use. Afterwards, you can drag-and-drop, or select one or multiple FASTQ or FASTA files to process them. Afterwards, you can get back your filtered (or selected) reads (or sequences). Below we also provide an example small human respiratory tract metagenome you can use.
            </p>

            <div>
              <p class="font-medium text-gray-900 mb-2">Example index and reads</p>
              <p class="text-sm">
                <a href="https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w61.idx" target="_blank" class="text-blue-600 hover:underline">Human read index</a>
              </p>
              <ul class="list-disc list-inside mt-1">
                <li>
                  <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/057/SRR31126957/SRR31126957_1.fastq.gz" class="text-blue-600 hover:underline">SRR31126957_1.fastq.gz</a>
                </li>
                <li>
                  <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR311/057/SRR31126957/SRR31126957_2.fastq.gz" class="text-blue-600 hover:underline">SRR31126957_2.fastq.gz</a>
                </li>
              </ul>
            </div>
          </TabsContent>

          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">Deplete mode</dt>
                <dd class="ml-4">Activated by default, this checkbox can be unticked to keep those reads that have matches with the index's, instead of removing them.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Absolute threshold</dt>
                <dd class="ml-4">This is the minimum number of k-mers that must match between a read and an index to select a read (and filter/keep it, depending on whether depleted mode is activated or not). Both this and the relative threshold must be satisfied to select a read.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Relative threshold</dt>
                <dd class="ml-4">This is the percentage (set as a real number between 0 and 1) of k-mers that must match between a read and an index to select a read (and filter/keep it, depending on whether depleted mode is activated or not). Both this and the absolute threshold must be satisfied to select a read.</dd>
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
  name: "HostDepletionHelpCollapsible",
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
