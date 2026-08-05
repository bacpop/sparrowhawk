<template>
  <div>
    <div class="flex items-center justify-between mb-2">
      <h2 class="text-base font-medium">Transmission clusters</h2>
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

    <div v-if="hasTimelineMatches" class="mt-6">
      <h3 class="text-sm font-medium text-gray-600 mb-2">Cluster timeline</h3>
      <TransmissionClusterTimeline
        :clusterRows="clusterTableRows"
        :metadataRows="metadataRows"
      />
    </div>
    <p v-else-if="metadataRows.length > 0" class="mt-6 text-sm text-gray-500">
      No metadata rows with matching sample IDs and dates were found for the cluster timeline.
    </p>
  </div>
</template>

<script lang="ts">
import { defineComponent, PropType } from "vue";
import { FileDown } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import TransmissionGraph from "@/components/TransmissionGraph.vue";
import TransmissionClusterTimeline from "@/components/TransmissionClusterTimeline.vue";
import { saveTextFile } from "@/platform/files";
import type { Dict, MetadataRow, TransmissionGraphData } from "@/types";

export default defineComponent({
  name: "TransmissionClusterResults",
  components: {
    Button,
    FileDown,
    TransmissionGraph,
    TransmissionClusterTimeline,
  },
  props: {
    clusterResults: {
      type: Object as PropType<Dict<number> | null>,
      required: false,
      default: null,
    },
    transmissionGraph: {
      type: Object as PropType<TransmissionGraphData | null>,
      required: false,
      default: null,
    },
    metadataRows: {
      type: Array as PropType<MetadataRow[]>,
      required: true,
    },
    enableLocationMatching: {
      type: Boolean,
      required: true,
    },
  },
  computed: {
    adjustedClusterResults(): Record<string, number> | null {
      if (!this.clusterResults || !this.enableLocationMatching || !this.metadataRows.length) {
        return this.clusterResults;
      }

      const keyToId = new Map<string, number>();
      let nextId = 1;
      const result: Record<string, number> = {};
      for (const [name, cluster] of Object.entries(this.clusterResults)) {
        const loc = this.sampleLocation(name) ?? '';
        const key = `${cluster}:${loc}`;
        if (!keyToId.has(key)) {
          keyToId.set(key, nextId++);
        }
        result[name] = keyToId.get(key)!;
      }
      return result;
    },
    clusterTableRows(): { name: string; cluster: number }[] {
      const results = this.adjustedClusterResults;
      if (!results) {
        return [];
      }
      return Object.entries(results)
        .map(([name, cluster]) => ({ name, cluster }))
        .sort((a, b) => a.cluster - b.cluster || a.name.localeCompare(b.name));
    },
    hasTimelineMatches(): boolean {
      return this.clusterTableRows.some((row) => {
        const [prefix] = row.name.split(".");
        return this.metadataRows.some((metadataRow) => (metadataRow.id === prefix || metadataRow.id === row.name) && Boolean(metadataRow.date));
      });
    },
    transmissionGraphNodes() {
      const nodes = this.transmissionGraph?.nodes ?? [];
      const adjusted = this.adjustedClusterResults;
      if (!adjusted) {
        return nodes;
      }
      return nodes.map((node) => ({
        ...node,
        cluster: adjusted[node.id] ?? node.cluster,
      }));
    },
    transmissionGraphLinks() {
      return this.transmissionGraph?.links ?? [];
    },
  },
  methods: {
    sampleLocation(name: string): string | null {
      const prefix = name.replace(/(?:\.[^.]+)+$/, '');
      const row = this.metadataRows.find((metadataRow) => metadataRow.id === prefix || metadataRow.id === name);
      return row?.location ?? null;
    },
    async downloadClustersTSV(): Promise<void> {
      const tsv = "Sample\tCluster\n" + this.clusterTableRows.map((row) => `${row.name}\t${row.cluster}`).join("\n");
      await saveTextFile(tsv, "transmission_clusters.tsv", "text/tab-separated-values;charset=utf-8");
    },
  },
});
</script>
