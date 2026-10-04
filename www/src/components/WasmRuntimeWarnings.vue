<!-- For showing a warning in case wasm64 binaries don't load -->

<template>
  <div v-if="warnings.length" class="mb-4 flex flex-col gap-2" role="status" aria-live="polite">
    <p v-for="notice in warnings" :key="notice.moduleId"
       class="rounded-md border border-amber-300 bg-amber-50 p-3 text-sm text-amber-900">
      {{ notice.label }}’s wasm64 binary failed to load. Using wasm32; larger datasets may encounter its memory limit.
    </p>
  </div>
</template>

<script lang="ts">
import { computed, defineComponent } from "vue";
import { useStore } from "vuex";
import type { RootState } from "@/store/state";

export default defineComponent({
  name: "WasmRuntimeWarnings",
  setup() {
    const store = useStore<RootState>();
    const warnings = computed(() => Object.values(store.state.wasmRuntimeStatus)
      .filter((notice) => notice.fallbackReason !== undefined));
    return { warnings };
  },
});
</script>
