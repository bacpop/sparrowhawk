<template>
  <div ref="plotEl" class="w-full" style="min-height: 300px;"></div>
</template>

<script lang="ts">
import { defineComponent, PropType } from "vue";
// @ts-ignore - plotly.js-dist doesn't have types
import Plotly from 'plotly.js-dist';

interface ClusterRow {
  name: string;
  cluster: number;
}

interface MetadataRow {
  id: string;
  date: string;
  location: string;
}

interface PreparedPoint {
  cluster: number;
  date: string;
  cluster_size: number;
  cases: number;
  text: string;
}

export default defineComponent({
  name: 'TransmissionClusterTimeline',

  props: {
    clusterRows: {
      type: Array as PropType<ClusterRow[]>,
      required: true,
    },
    metadataRows: {
      type: Array as PropType<MetadataRow[]>,
      required: true,
    },
  },

  computed: {
    preparedData(): PreparedPoint[] {
      // Join cluster rows with metadata by matching sample name (strip file extensions)
      const joined: { cluster: number; date: string }[] = [];
      for (const row of this.clusterRows) {
        const prefix = row.name.replace(/(?:\.[^.]+)+$/, '');
        const meta = this.metadataRows.find(m => m.id === prefix || m.id === row.name);
        if (meta?.date) {
          joined.push({ cluster: row.cluster, date: meta.date });
        }
      }
      if (joined.length === 0) return [];

      // cluster_size: count per cluster
      const clusterSize = new Map<number, number>();
      for (const { cluster } of joined) {
        clusterSize.set(cluster, (clusterSize.get(cluster) ?? 0) + 1);
      }

      // cases: count per (cluster, date)
      const casesMap = new Map<string, number>();
      for (const { cluster, date } of joined) {
        const key = `${cluster}::${date}`;
        casesMap.set(key, (casesMap.get(key) ?? 0) + 1);
      }

      // Deduplicate to one row per (cluster, date)
      const seen = new Set<string>();
      const result: PreparedPoint[] = [];
      for (const { cluster, date } of joined) {
        const key = `${cluster}::${date}`;
        if (seen.has(key)) continue;
        seen.add(key);
        const cluster_size = clusterSize.get(cluster) ?? 0;
        const cases = casesMap.get(key) ?? 0;
        const text = `Cluster ${cluster}\nCluster size: ${cluster_size}\nDate: ${date}\nCases: ${cases}`;
        result.push({ cluster, date, cluster_size, cases, text });
      }

      return result.sort((a, b) => a.cluster - b.cluster || a.date.localeCompare(b.date));
    },
  },

  methods: {
    renderChart(): void {
      const data = this.preparedData;
      if (!data.length || !this.$refs.plotEl) return;

      const clusters = [...new Set(data.map(d => d.cluster))].sort((a, b) => a - b);
      const traces: object[] = [];

      for (const cl of clusters) {
        const points = data.filter(d => d.cluster === cl).sort((a, b) => a.date.localeCompare(b.date));
        const xs = points.map(p => p.date);
        const ys = points.map(() => String(cl));

        // Line trace connecting points within the cluster
        traces.push({
          x: xs,
          y: ys,
          type: 'scatter',
          mode: 'lines',
          line: { color: 'grey', width: 1.8 },
          opacity: 0.75,
          showlegend: false,
          hoverinfo: 'skip',
        });

        // Marker trace sized by cases
        traces.push({
          x: xs,
          y: ys,
          type: 'scatter',
          mode: 'markers',
          marker: {
            size: points.map(p => p.cases * 5),
            sizemode: 'diameter',
            sizeref: 1,
            showscale: false,
          },
          text: points.map(p => p.text),
          hoverinfo: 'text',
          name: `Cluster ${cl}`,
          showlegend: false,
        });
      }

      const layout = {
        xaxis: { title: 'Isolation date', type: 'date', tickformat: '%Y-%m', tickangle: 45 },
        yaxis: { title: 'Cluster', type: 'category' },
        showlegend: false,
        font: { family: 'IBM Plex Sans' },
        autosize: true,
        margin: { t: 30, r: 20 },
      };

      Plotly.newPlot(this.$refs.plotEl as HTMLElement, traces, layout as Partial<Plotly.Layout>);
    },
  },

  watch: {
    preparedData: {
      handler(): void {
        this.$nextTick(() => this.renderChart());
      },
      deep: true,
    },
  },
});
</script>
