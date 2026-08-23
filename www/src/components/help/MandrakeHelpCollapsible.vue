<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Mandrake</span>
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
          <TabsContent value="overview" class="space-y-3">
            <p>
              Upload an alignment, a Roary-style accessory table, or a paired
              <code>.skm</code>/<code>.skd</code> sketch database. Mandrake builds sparse
              distances and calculates a two-dimensional stochastic cluster embedding locally.
            </p>
            <p>
              The embedding appears progressively. You can optionally provide a two-column
              sample-name/label TSV to colour the plot, run HDBSCAN on the final coordinates,
              and download the coordinates, names, or cluster labels.
            </p>
          </TabsContent>
          <TabsContent value="parameters">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">Input formats</dt>
                <dd class="ml-4">FASTA/FASTQ inputs are treated as aligned sequences; Rtab/TSV inputs are accessory profiles. Sketch databases require both paired files.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Sparsification</dt>
                <dd class="ml-4">Keep a fixed number of nearest neighbours or retain distances below a strict normalized threshold.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">HDBSCAN</dt>
                <dd class="ml-4">Optional fixed-preset clustering of the completed two-dimensional embedding; clustering failures are reported without discarding the embedding.</dd>
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
import { ChevronDown, HelpCircle } from "@lucide/vue";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@/components/ui/collapsible";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

export default defineComponent({
  name: "MandrakeHelpCollapsible",
  components: {
    ChevronDown,
    HelpCircle,
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
