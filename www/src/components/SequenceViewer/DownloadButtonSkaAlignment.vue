<template>
  <div class="flex flex-col gap-3">
    <div class="flex flex-wrap gap-6 items-start">
      <div class="flex flex-col gap-2">
        <span>Download alignment</span>
        <div v-if="hasAlignment" class="flex flex-wrap items-end gap-2">
          <Button v-if="supportsStreaming" :disabled="isBusy || isDownloading || alignmentExportComplete" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadALN">
            <Loader2 v-if="isDownloading" class="w-4 h-4 mr-2 animate-spin" />
            <Download v-else class="w-4 h-4 mr-2" />
            .aln.gz
          </Button>
          <template v-else>
            <Button :disabled="isBusy || isDownloading || (alignmentExportComplete && !currentPartBlob)" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadCurrentPart">
              <Loader2 v-if="isDownloading" class="w-4 h-4 mr-2 animate-spin" />
              <Download v-else class="w-4 h-4 mr-2" />
              {{ totalParts > 1 ? `.aln.gz (${partLabel(currentPartIndex)})` : '.aln.gz' }}
            </Button>
            <Button v-if="totalParts > 1 && hasNextPart" :disabled="isBusy || isDownloading || !currentPartBlob" class="max-w-fit cursor-pointer" variant="outline" size="sm" @click="downloadNextPart">
              <Download class="w-4 h-4 mr-2" />
              .aln.gz ({{ partLabel(currentPartIndex + 1) }})
            </Button>
          </template>
        </div>
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
    <p v-if="!supportsStreaming && totalParts > 1" class="text-xs text-gray-500">
      Download each part in order. The current part stays available until you request the next one.
    </p>
    <p v-if="downloadError" class="text-sm text-red-600" role="alert">{{ downloadError }}</p>
  </div>
</template>

<script lang="ts">
import { defineComponent } from "vue";
import { useStore } from "vuex";
import { Download, Loader2 } from "@lucide/vue";
import { Button } from "@/components/ui/button";
import { openBinaryFileSink, saveBinaryFile, saveBlobFile, saveTextFile, supportsStreamingFileSave } from "@/platform/files";
import type { Alignment } from "@/types";
import type { RootState } from "@/store/state";

const ALIGNMENT_EXPORT_SAMPLES_PER_PART = 256;

export default defineComponent({
  name: 'DownloadButtonSkaAlignment',
  components: { Button, Download, Loader2 },
  setup() { return { store: useStore<RootState>() }; },
  data() {
    return {
      isDownloading: false,
      saveError: "",
      currentPartBlob: null as Blob | null,
      currentPartIndex: 1,
      exportStarted: false,
    };
  },
  computed: {
    result(): Alignment | undefined { return this.store.state.allResults_ska.alignResults[0]; },
    hasAlignment(): boolean { return this.result?.alignmentAvailable === true; },
    alignmentExportComplete(): boolean { return this.result?.alignmentExportComplete === true; },
    totalSamples(): number { return this.result?.names?.length ?? 0; },
    totalParts(): number { return Math.ceil(this.totalSamples / ALIGNMENT_EXPORT_SAMPLES_PER_PART); },
    supportsStreaming(): boolean { return supportsStreamingFileSave(); },
    hasNextPart(): boolean { return this.currentPartIndex < this.totalParts; },
    hasDistancesGzip(): boolean { return (this.result?.distances_csv_gzip?.byteLength ?? 0) > 0; },
    downloadError(): string { return this.saveError || this.result?.alignmentDownloadError || ""; },
    isBusy(): boolean {
      const p = this.store.state.processingState;
      return p.isAligning || p.isObtainingAlignment || p.isClustering ||
        p.isTransmissionStandaloneClustering || p.isIndexingRef || p.isMapping;
    },
  },
  watch: {
    result(next: Alignment | undefined, previous: Alignment | undefined): void {
      if (next?.runId === previous?.runId) return;
      this.currentPartBlob = null;
      this.currentPartIndex = 1;
      this.exportStarted = false;
      this.saveError = "";
    },
  },
  methods: {
    async downloadALN(): Promise<void> {
      if (this.isBusy || this.isDownloading) return;
      this.isDownloading = true;
      this.saveError = "";
      let sink: Awaited<ReturnType<typeof openBinaryFileSink>> = null;
      try {
        sink = await openBinaryFileSink("alignment.aln.gz", "application/gzip");
        if (!sink) return;
        await this.store.dispatch("streamAlignmentExport", { sink });
      } catch (error) {
        if (error instanceof DOMException && error.name === "AbortError") return;
        this.saveError = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment download failed", error);
      } finally { this.isDownloading = false; }
    },
    ordinal(value: number): string {
      const suffix = value % 100 >= 11 && value % 100 <= 13 ? "th" : ({ 1: "st", 2: "nd", 3: "rd" } as Record<number, string>)[value % 10] ?? "th";
      return `${value}${suffix}`;
    },
    partLabel(index: number): string {
      return `${this.ordinal(index)} part${index === this.totalParts ? ", last" : ""}`;
    },
    async generatePart(start: boolean): Promise<Blob> {
      const blob: Blob = await this.store.dispatch("getAlignmentExportPart", { samplesPerPart: ALIGNMENT_EXPORT_SAMPLES_PER_PART, start });
      this.exportStarted = true;
      return blob;
    },
    saveCurrentPart(): void {
      if (!this.currentPartBlob) return;
      if (this.totalSamples <= ALIGNMENT_EXPORT_SAMPLES_PER_PART) {
        saveBlobFile(this.currentPartBlob, "alignment.aln.gz");
        return;
      }
      const part = String(this.currentPartIndex).padStart(3, "0");
      saveBlobFile(this.currentPartBlob, `alignment.part-${part}.aln.gz`);
    },
    async downloadCurrentPart(): Promise<void> {
      if (this.isBusy || this.isDownloading) return;
      this.isDownloading = true;
      this.saveError = "";
      try {
        if (!this.currentPartBlob) this.currentPartBlob = await this.generatePart(!this.exportStarted);
        this.saveCurrentPart();
      } catch (error) {
        this.currentPartBlob = null;
        this.currentPartIndex = 1;
        this.exportStarted = false;
        this.saveError = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment part download failed", error);
      } finally { this.isDownloading = false; }
    },
    async downloadNextPart(): Promise<void> {
      if (this.isBusy || this.isDownloading || !this.currentPartBlob || !this.hasNextPart) return;
      this.isDownloading = true;
      this.saveError = "";
      this.currentPartBlob = null;
      this.currentPartIndex += 1;
      try {
        this.currentPartBlob = await this.generatePart(false);
        this.saveCurrentPart();
      } catch (error) {
        // The worker resets its cursors on cancellation; restart from part one on retry.
        this.currentPartBlob = null;
        this.currentPartIndex = 1;
        this.exportStarted = false;
        this.saveError = error instanceof Error ? error.message : String(error);
        console.error("[ska-align] Alignment part download failed", error);
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
