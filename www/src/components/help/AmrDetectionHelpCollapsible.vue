<template>
  <div class="mb-4 rounded-md border border-gray-200 bg-gray-50 p-4 text-sm text-gray-700">
    AMR detection uses a browser-side Rust/WebAssembly index built from ResFinderDB.
    Upload a prebuilt AMR index and one or more FASTA files to detect acquired resistance genes locally.
  </div>
</template>

<script lang="ts">
import { defineComponent } from "vue";

export default defineComponent({
  name: "AmrDetectionHelpCollapsible",
});
</script>
