<template>
  <div class="flex flex-col gap-3">
    <div class="flex flex-wrap gap-6 items-start">
      <div class="flex flex-col gap-2">
        <span>Download alignment</span>
        <Button v-if="hasAlignment" :disabled="isBusy || isDownloading" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadALN">
          <Loader2 v-if="isDownloading" class="w-4 h-4 mr-2 animate-spin" />
          <Download v-else class="w-4 h-4 mr-2" />
          .aln.gz
        </Button>
      </div>
      <div class="flex flex-col gap-2">
        <span>Download tree</span>
        <Button v-if="hasAlignment" :disabled="isBusy" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadTREE">
          <Download class="w-4 h-4 mr-2" /> .tree
        </Button>
      </div>
      <div class="flex flex-col gap-2">
        <span>Download distances</span>
        <Button v-if="hasDistancesGzip" :disabled="isBusy" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadCSV">
          <Download class="w-4 h-4 mr-2" /> .csv.gz
        </Button>
      </div>
    </div>
    <p v-if="hasAlignment && !result?.alignmentFrozen" class="text-xs text-gray-500">
      The first alignment download closes this dataset to new samples.
    </p>
    <p v-if="downloadError" class="text-sm text-red-600" role="alert">{{ downloadError }}</p>
  </div>
</template>

<script lang="ts">
import { defineComponent } from "vue";
import { useStore } from "vuex";
import { Download, Loader2 } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import { saveBinaryFile, saveTextFile } from "@/platform/files";
import type { Alignment } from "@/types";
import type { RootState } from "@/store/state";

export default defineComponent({
  name: 'DownloadButtonSkaAlignment',
  components: { Button, Download, Loader2 },
  setup() { return { store: useStore<RootState>() }; },
  data() { return { isDownloading: false, saveError: "" }; },
  computed: {
    result(): Alignment | undefined { return this.store.state.allResults_ska.alignResults[0]; },
    hasAlignment(): boolean { return this.result?.alignmentAvailable === true; },
    hasDistancesGzip(): boolean { return (this.result?.distances_csv_gzip?.byteLength ?? 0) > 0; },
    downloadError(): string { return this.saveError || this.result?.alignmentDownloadError || ""; },
    isBusy(): boolean {
      const p = this.store.state.processingState;
      return p.isAligning || p.isObtainingAlignment || p.isClustering ||
        p.isTransmissionStandaloneClustering || p.isIndexingRef || p.isMapping;
    },
  },
  methods: {
    async downloadALN(): Promise<void> {
      if (this.isBusy || this.isDownloading) return;
      this.isDownloading = true;
      this.saveError = "";
      try {
        const bytes: Uint8Array = await this.store.dispatch("getAlignmentDownload");
        await saveBinaryFile(bytes, "alignment.aln.gz", "application/gzip");
      } catch (error) {
        this.saveError = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment download failed", error);
      } finally { this.isDownloading = false; }
    },
    async downloadTREE(): Promise<void> {
      if (this.result?.newick) await saveTextFile(this.result.newick, "alignmenttree.tree", "text/x-newick");
    },
    async downloadCSV(): Promise<void> {
      const bytes = this.result?.distances_csv_gzip;
      if (bytes?.byteLength) await saveBinaryFile(bytes, "distances.csv.gz", "application/gzip");
    },
  },
});
</script>
