<template>
  <Collapsible v-model:open="isOpen" class="border border-gray-200 rounded-lg bg-white mb-4">
    <CollapsibleTrigger class="flex items-center justify-between w-full p-3 hover:bg-gray-50 rounded-lg">
      <div class="flex items-center gap-2">
        <HelpCircle class="w-4 h-4 text-gray-500" />
        <span class="font-medium text-sm">How to use Transmission</span>
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
          <TabsTrigger value="inputs" class="text-xs">Inputs</TabsTrigger>
        </TabsList>

        <div class="mt-3 text-sm text-gray-600 max-h-64 overflow-y-auto">
          <TabsContent value="overview" class="space-y-4">
            <p>
              Use this tab to run transmission clustering from an alignment that has already been generated elsewhere.
              The clustering uses single-linkage components under the selected SNP threshold.
            </p>
            <p>
              SNP distances are computed from aligned columns. Positions with gaps, Ns, or ambiguity codes in either
              sample are skipped. U is treated as T.
            </p>
          </TabsContent>

          <TabsContent value="inputs">
            <dl class="space-y-3">
              <div>
                <dt class="font-medium text-gray-900">Alignment</dt>
                <dd class="ml-4">Upload an aligned FASTA file. Files ending in .aln are accepted when their content is FASTA-formatted.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">SNP threshold</dt>
                <dd class="ml-4">Maximum pairwise SNP distance for two samples to be connected in the transmission graph.</dd>
              </div>
              <div>
                <dt class="font-medium text-gray-900">Metadata</dt>
                <dd class="ml-4">Optionally upload a CSV with ID, date, and location columns to show timelines and split clusters by location.</dd>
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
import { HelpCircle, ChevronDown } from "lucide-vue-next";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@/components/ui/collapsible";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

export default defineComponent({
  name: "TransmissionHelpCollapsible",
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
