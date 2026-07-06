<template>
  <div class="flex flex-col gap-2">
    <Button @click="downloadR1" class="mt-4" variant="outline">
      <Download class="w-4 h-4 mr-2" />
      Download {{ isPaired ? 'R1' : 'filtered reads' }} (.fastq.gz)
    </Button>
    <Button v-if="isPaired" @click="downloadR2" variant="outline">
      <Download class="w-4 h-4 mr-2" />
      Download R2 (.fastq.gz)
    </Button>
  </div>
</template>
<script lang="ts">
import { defineComponent, computed } from "vue";
import { useStore } from "vuex";
import { Button } from "@/components/ui/button";
import { Download } from "lucide-vue-next";
import {saveBinaryFile} from "@/platform/files";
export default defineComponent({
  name: "DownloadButtonHostDepletion",
  components: { Button, Download },
  setup() {
    const store = useStore();

    const isPaired = computed(() => store.state.allResults_deacon.outputGzip2 !== null);

    async function downloadFile(content: Uint8Array<ArrayBuffer>, baseName: string): Promise<void> {
      await saveBinaryFile(content, baseName, "application/gzip");
    }

    async function downloadR1(): Promise<void> {
      const content: Uint8Array | null = store.state.allResults_deacon.outputGzip;
      if (!content) return;
      const base = (store.state.allResults_deacon.readsFileName ?? "filtered")
                     .replace(/\.(fastq|fq)(\.gz)?$/, '');
      await downloadFile(new Uint8Array(content), base + "_filtered.fastq.gz");
    }

    async function downloadR2(): Promise<void> {
      const content: Uint8Array | null = store.state.allResults_deacon.outputGzip2;
      if (!content) return;
      const base = (store.state.allResults_deacon.readsFileName2 ?? "filtered_R2")
                     .replace(/\.(fastq|fq)(\.gz)?$/, '');
      await downloadFile(new Uint8Array(content), base + "_filtered.fastq.gz");
    }

    return { isPaired, downloadR1, downloadR2 };
  }
});
</script>
