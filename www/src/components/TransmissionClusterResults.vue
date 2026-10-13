<template>
  <div>
    <div class="flex items-center justify-between mb-2">
      <h2 class="text-base font-medium">Transmission clusters</h2>
    </div>

    <p v-if="downloadError" class="mb-2 text-sm text-red-600" role="alert">
      {{ downloadError }}
    </p>

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
            <td class="px-3 py-2" :class="{ 'text-gray-400': row.singleton }">
              {{ row.cluster }}<span v-if="row.singleton"> (singleton)</span>
            </td>
          </tr>
        </tbody>
      </table>
    </div>

    <h2 class="mt-6 mb-2 text-base font-medium">Network graph</h2>
    <div v-if="transmissionGraphNodes.length" class="h-[500px]">
      <TransmissionGraph
        :nodes="transmissionGraphNodes"
        :links="transmissionGraphLinks"
      >
        <template v-if="singletonCount > 0" #toolbar-note>
          {{ singletonCount }} sample(s) without links at this SNP threshold are omitted from the graph and timeline;
          they are shown in grey in the tree.
        </template>
      </TransmissionGraph>
    </div>
    <p v-else class="text-sm text-gray-500">No samples are linked at the selected SNP threshold.</p>

    <div v-if="hasTimelineMatches" class="mt-6">
      <h2 class="text-base font-medium mb-2">Cluster timeline</h2>
      <TransmissionClusterTimeline
        :clusterRows="clusterTableRows"
        :metadataRows="metadataRows"
      />
    </div>
    <p v-else-if="metadataRows.length > 0" class="mt-6 text-sm text-gray-500">
      No linked samples with matching metadata IDs and dates were found for the cluster timeline.
    </p>
  </div>
</template>

<script lang="ts">
import { defineComponent, PropType } from "vue";
import TransmissionGraph from "@/components/TransmissionGraph.vue";
import TransmissionClusterTimeline from "@/components/TransmissionClusterTimeline.vue";
import { buildSampleMetadataIndex, adjustClusterAssignments, singletonSampleNames } from "@/utils/transmissionMetadata";
import type { SampleMetadataIndex } from "@/utils/transmissionMetadata";
import { saveTextFile } from "@/platform/files";
import type { Dict, MetadataRow, TransmissionGraphData } from "@/types";

export default defineComponent({
  name: "TransmissionClusterResults",
  components: {
    TransmissionGraph,
    TransmissionClusterTimeline,
  },
  emits: ["download-state"],
  data() { return { isDownloading: false, downloadError: "" }; },
  props: {
    downloadDisabled: {
      type: Boolean,
      default: false,
    },
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
    metadataIndex(): SampleMetadataIndex {
      return buildSampleMetadataIndex(Object.keys(this.clusterResults ?? {}), this.metadataRows);
    },
    adjustedClusterResults(): Record<string, number> | null {
      return adjustClusterAssignments(this.clusterResults, this.metadataIndex, this.enableLocationMatching);
    },
    singletonNames(): Set<string> {
      return singletonSampleNames(this.clusterResults);
    },
    singletonCount(): number {
      return this.singletonNames.size;
    },
    clusterTableRows(): { name: string; cluster: number; singleton: boolean }[] {
      const results = this.adjustedClusterResults;
      if (!results) {
        return [];
      }
      return Object.entries(results)
        .map(([name, cluster]) => ({ name, cluster, singleton: this.singletonNames.has(name) }))
        .sort((a, b) => a.cluster - b.cluster || a.name.localeCompare(b.name));
    },
    hasTimelineMatches(): boolean {
      return this.clusterTableRows.some(row =>
        !row.singleton && Boolean(this.metadataIndex.rowsBySample.get(row.name)?.date));
    },
    transmissionGraphNodes() {
      const adjusted = this.adjustedClusterResults;
      return (this.transmissionGraph?.nodes ?? [])
        .filter((node) => !this.singletonNames.has(node.id))
        .map((node) => ({
          ...node,
          cluster: adjusted?.[node.id] ?? node.cluster,
        }));
    },
    transmissionGraphLinks() {
      return this.transmissionGraph?.links ?? [];
    },
  },
  methods: {
    async downloadClustersTSV(): Promise<void> {
      if (this.downloadDisabled || this.isDownloading || !this.clusterTableRows.length) return;

      this.isDownloading = true;
      this.downloadError = "";
      this.$emit("download-state", true);
      try {
        const tsv = "Sample\tCluster\n" + this.clusterTableRows
          .map((row) => `${row.name}\t${row.cluster}`).join("\n");
        await saveTextFile(tsv, "transmission_clusters.tsv", "text/tab-separated-values;charset=utf-8");
      } catch (error) {
        this.downloadError = error instanceof Error ? error.message : String(error);
        console.error("[transmission] Cluster download failed", error);
      } finally {
        this.isDownloading = false;
        this.$emit("download-state", false);
      }
    },
  },
});
</script>
