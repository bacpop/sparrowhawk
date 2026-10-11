<template>
  <TaxoniumTree v-if="tree" :tree="tree" :cluster-results="clusterResults"
    :metadata-rows="metadataRows" :enable-location-matching="enableLocationMatching" />
</template>

<script lang="ts">
import { computed, defineComponent, PropType } from 'vue';
import { useStore } from 'vuex';
import TaxoniumTree from '@/components/TaxoniumTree.vue';
import type { RootState } from '@/store/state';
import type { MetadataRow, TreeResult } from '@/types';

export default defineComponent({
  name: 'ResultsDisplayAlignment',
  components: { TaxoniumTree },
  props: {
    metadataRows: { type: Array as PropType<MetadataRow[]>, default: () => [] },
    enableLocationMatching: { type: Boolean, default: false },
  },
  setup() {
    const store = useStore<RootState>();
    const tree = computed<TreeResult | null>(() => {
      const alignment = store.state.allResults_ska.alignResults[0];
      if (!alignment?.aligned || alignment.invalidated) return null;
      return { names: alignment.names ?? [], treeTipNames: alignment.treeTipNames, newick: alignment.newick ?? '' };
    });
    const clusterResults = computed(() => store.state.allResults_ska.clusterResults);
    return { tree, clusterResults };
  },
});
</script>
