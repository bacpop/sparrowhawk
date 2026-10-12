<template>
  <div v-if="tree" class="min-w-0 w-full">
    <p v-if="tree.names.length < 3" class="text-sm text-gray-500" role="status">
      At least three samples are required to visualise a tree.
    </p>
    <template v-else-if="tree.newick">
      <p v-if="!hasTipMapping" class="mb-2 text-sm text-amber-800" role="status">
        Re-run analysis to add tree metadata.
      </p>
      <p v-else-if="metadataIndex.ambiguousSamples.length" class="mb-2 text-sm text-amber-800" role="status">
        Metadata IDs are ambiguous for {{ metadataIndex.ambiguousSamples.length }} sample(s); their date and location fields are blank.
      </p>
      <div ref="treeContainer" class="w-full min-w-0" style="height: 550px;"></div>
    </template>
    <p v-else class="text-sm text-gray-500" role="status">No tree is available.</p>
  </div>
</template>



<script lang="ts">
    import { computed, defineComponent, onBeforeUnmount, PropType, ref, watch } from 'vue';
    import * as ReactDOM from 'react-dom/client';
    import * as React from 'react';
    import Taxonium from 'taxonium-component';
    import type { Dict, MetadataRow, TreeResult } from '@/types';
    import { adjustClusterAssignments, buildSampleMetadataIndex } from '@/utils/transmissionMetadata';
    import { buildTaxoniumClusterColorMapping } from '@/utils/transmissionClusterColors';


    export default defineComponent({
      name: 'TaxoniumTree',

      props: {
        tree: { type: Object as PropType<TreeResult | null>, default: null },
        clusterResults: { type: Object as PropType<Dict<number> | null>, default: null },
        metadataRows: { type: Array as PropType<MetadataRow[]>, default: () => [] },
        enableLocationMatching: { type: Boolean, default: false },
      },

      setup(props) {
        const treeContainer = ref<HTMLDivElement | null>(null);

        // These will change whenever the arguments change, it's equivalent to put them when you define a website in the "computed" section
        const metadataIndex = computed(() => buildSampleMetadataIndex(props.tree?.names ?? [], props.metadataRows));
        const clusters      = computed(() => adjustClusterAssignments(props.clusterResults, metadataIndex.value, props.enableLocationMatching));
        const hasTipMapping = computed(() => {
          const tree = props.tree;
          const tips = tree?.treeTipNames;
          return Boolean(tree && tips && tips.length === tree.names.length &&
            new Set(tree.names).size === tree.names.length && new Set(tips).size === tips.length &&
            tips.every(tip => /^[A-Za-z0-9._-]+$/.test(tip)));
        });

        const metadataTsv = computed(() => {
          if (!props.tree || !hasTipMapping.value) return '';
          const clean = (value: string) => value.replace(/[\t\r\n]/g, ' ');
          const rows = props.tree.names.map((name, i) => {
            const metadata = metadataIndex.value.rowsBySample.get(name);
            const cluster = clusters.value?.[name];
            return [props.tree!.treeTipNames![i], cluster == null ? '' : `Cluster ${cluster}`,
              name, metadata?.date ?? '', metadata?.location ?? ''].map(clean).join('\t');
          });
          return ['Sample\tCluster\tOriginal sample\tDate\tLocation', ...rows].join('\n');
        });

        const input = computed(() => {
          const tree = props.tree;
          if (!tree?.newick || tree.names.length < 3) return null;
          const sourceData = {
            status: 'loaded', filename: 'alignment.nwk', filetype: 'nwk', data: tree.newick,
            ...(metadataTsv.value ? { metadata: {
              status: 'loaded', filename: 'alignment_metadata.tsv', filetype: 'meta_tsv', data: metadataTsv.value,
            } } : {}),
          };
          const defaultColorByField = hasTipMapping.value && tree.names.some(name => clusters.value?.[name] != null)
            ? 'meta_Cluster' : 'None';
          const colorMapping = buildTaxoniumClusterColorMapping(Object.values(clusters.value ?? {}));
          return { sourceData, configDict: { defaultColorByField, colorMapping },
            signature: JSON.stringify([tree.newick, metadataTsv.value, defaultColorByField, colorMapping]) };
        });

        // This is Taxonium stuff, that uses React
        let reactRoot: ReturnType<typeof ReactDOM.createRoot> | null = null;
        let rootContainer: HTMLDivElement | null = null;
        let resizeObserver: ResizeObserver | null = null;
        let frame: number | null = null;
        let disposed = false;
        let signature: string | null = null;
        let revision = 0;

        function unmountReact(): void {
          reactRoot?.unmount();
          reactRoot = null;
          rootContainer = null;
          signature = null;
        }

        function mountOrUpdate(): void {
          if (disposed) return;
          const container = treeContainer.value;
          const currentInput = input.value;
          if (rootContainer !== container || !currentInput) unmountReact();
          if (!container || !currentInput || container.clientWidth <= 0 || container.clientHeight <= 0) return;
          if (signature === currentInput.signature && reactRoot) return;
          if (!reactRoot) {
            reactRoot = ReactDOM.createRoot(container);
            rootContainer = container;
          }
          signature = currentInput.signature;
          // Taxonium loads its configuration on mount: metadata changes require a new component.
          reactRoot.render(React.createElement(Taxonium as React.ElementType, {
            sourceData: currentInput.sourceData, configDict: currentInput.configDict, key: ++revision,
          }));
        }

        function scheduleMount(): void {
          if (disposed || frame !== null) return;
          frame = requestAnimationFrame(() => { frame = null; mountOrUpdate(); });
        }

        watch(treeContainer, container => {
          resizeObserver?.disconnect();
          resizeObserver = null;
          if (rootContainer !== container) unmountReact();
          if (container) {
            resizeObserver = new ResizeObserver(scheduleMount);
            resizeObserver.observe(container);
          }
          mountOrUpdate();
        }, { flush: 'post' });
        watch(input, mountOrUpdate, { flush: 'post', immediate: true });

        onBeforeUnmount(() => {
          disposed = true;
          resizeObserver?.disconnect();
          if (frame !== null) cancelAnimationFrame(frame);
          unmountReact();
        });

        return { treeContainer, hasTipMapping, metadataIndex };
      },
    });
</script>
