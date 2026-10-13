<template>
  <div class="flex flex-col w-full h-full border border-gray-200 rounded-md bg-white overflow-hidden">
    <div class="shrink-0 flex items-center justify-between gap-3 border-b border-gray-200 p-2">
      <div class="min-w-0 text-sm text-gray-500"><slot name="toolbar-note" /></div>
      <Button type="button" variant="outline" size="sm" class="shrink-0" :disabled="!canFit" @click="fitGraph">
        Fit graph
      </Button>
    </div>
    <div ref="viewportRef" class="relative min-h-0 flex-1">
      <svg ref="svgRef" class="block w-full h-full" aria-label="Transmission cluster graph" />
      <p v-if="nodes.length === 0" class="absolute inset-0 flex items-center justify-center text-sm text-gray-500" role="status">
        No samples to display
      </p>
      <div v-if="clusterIds.length > 0"
           class="absolute top-3 right-3 bg-white bg-opacity-90 border border-gray-200 rounded-md p-2 text-xs space-y-1 max-h-48 overflow-y-auto">
        <div v-for="clusterId in clusterIds" :key="clusterId" class="flex items-center gap-1.5">
          <span class="inline-block w-3 h-3 rounded-full flex-shrink-0" :style="{ backgroundColor: clusterColor(clusterId) }" />
          <span>Cluster {{ clusterId }}</span>
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, ref, watch, onMounted, onBeforeUnmount, computed, PropType } from "vue";
import * as d3 from "d3";
import { Button } from "@/components/ui/button";
import type { TransmissionGraphNode, TransmissionGraphLink } from "@/types";
import { transmissionClusterColor } from "@/utils/transmissionClusterColors";

interface LayoutNode extends TransmissionGraphNode, d3.SimulationNodeDatum {}

interface LayoutLink extends d3.SimulationLinkDatum<LayoutNode> {
  source: LayoutNode;
  target: LayoutNode;
  snp_distance: number;
}

interface GraphBounds {
  x: number;
  y: number;
  width: number;
  height: number;
}

export default defineComponent({
  name: "TransmissionGraph",
  components: { Button },
  props: {
    nodes: {
      type: Array as PropType<TransmissionGraphNode[]>,
      required: true,
    },
    links: {
      type: Array as PropType<TransmissionGraphLink[]>,
      required: true,
    },
  },
  setup(props) {
    const viewportRef = ref<HTMLDivElement | null>(null);
    const svgRef = ref<SVGSVGElement | null>(null);
    const canFit = ref(false);
    let simulation: d3.Simulation<LayoutNode, LayoutLink> | null = null;
    let zoomBehavior: d3.ZoomBehavior<SVGSVGElement, unknown> | null = null;
    let contentGroup: d3.Selection<SVGGElement, unknown, null, undefined> | null = null;
    let layoutNodes: LayoutNode[] = [];
    let cameraMode: "fit" | "manual" = "fit";
    let viewportWidth = 0;
    let viewportHeight = 0;
    let viewportVisible = false;
    let needsBuild = true;
    let disposed = false;
    let resizeObserver: ResizeObserver | null = null;
    let resizeFrame: number | null = null;

    const clusterIds = computed(() => {
      const ids = new Set(props.nodes.map(n => n.cluster));
      return Array.from(ids).sort((a, b) => a - b);
    });

    const clusterColor = transmissionClusterColor;

    function hasFinitePositions(): boolean {
      return layoutNodes.some(node => Number.isFinite(node.x) && Number.isFinite(node.y));
    }

    function getGraphBounds(): GraphBounds | null {
      if (!hasFinitePositions()) return null;
      try {
        const bounds = contentGroup?.node()?.getBBox();
        if (bounds && Number.isFinite(bounds.x) && Number.isFinite(bounds.y)
            && Number.isFinite(bounds.width) && Number.isFinite(bounds.height)
            && bounds.width > 0 && bounds.height > 0) {
          return { x: bounds.x, y: bounds.y, width: bounds.width, height: bounds.height };
        }
      } catch {
        // Hidden SVGs may not expose a bounding box; node coordinates still allow fitting.
      }

      let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
      for (const node of layoutNodes) {
        const { x, y } = node;
        if (x === undefined || y === undefined || !Number.isFinite(x) || !Number.isFinite(y)) continue;
        minX = Math.min(minX, x - 12);
        minY = Math.min(minY, y - 12);
        maxX = Math.max(maxX, x + 12);
        maxY = Math.max(maxY, y + 12);
      }
      return Number.isFinite(minX)
        ? { x: minX, y: minY, width: maxX - minX, height: maxY - minY }
        : null;
    }

    function applyTransform(transform: d3.ZoomTransform): void {
      if (disposed || !svgRef.value || !zoomBehavior) return;
      if (!Number.isFinite(transform.x) || !Number.isFinite(transform.y)
          || !Number.isFinite(transform.k) || transform.k <= 0) return;
      d3.select(svgRef.value).call(zoomBehavior.transform, transform);
    }

    function fitGraph(): void {
      if (disposed || !viewportVisible || !zoomBehavior) return;
      const bounds = getGraphBounds();
      if (!bounds) return;
      const paddingX = Math.min(50, viewportWidth / 4);
      const paddingY = Math.min(50, viewportHeight / 4);
      const scale = Math.min(
        (viewportWidth - 2 * paddingX) / Math.max(bounds.width, 1),
        (viewportHeight - 2 * paddingY) / Math.max(bounds.height, 1),
        3,
      );
      if (!Number.isFinite(scale) || scale <= 0) return;
      cameraMode = "fit";
      // Fit and interactive zoom must share limits, including fits below 0.2.
      zoomBehavior.scaleExtent([Math.min(0.2, scale / 2), 5]);
      applyTransform(d3.zoomIdentity.translate(
        viewportWidth / 2 - scale * (bounds.x + bounds.width / 2),
        viewportHeight / 2 - scale * (bounds.y + bounds.height / 2),
      ).scale(scale));
    }

    function disposeGraph(): void {
      simulation?.stop().on("tick", null).on("end", null);
      simulation = null;
      if (svgRef.value) {
        const svg = d3.select(svgRef.value);
        svg.interrupt().on(".zoom", null);
        svg.selectAll("*").remove();
      }
      contentGroup = null;
      zoomBehavior = null;
      layoutNodes = [];
      canFit.value = false;
    }

    function buildGraph(): void {
      if (disposed || !svgRef.value || !viewportVisible) return;
      needsBuild = false;
      cameraMode = "fit";
      if (props.nodes.length === 0) return;

      const svg = d3.select(svgRef.value);
      // D3 owns these copies and seeds distributed positions when x/y are absent.
      layoutNodes = props.nodes.map(node => ({ ...node, fx: null, fy: null }));
      const nodeById = new Map(layoutNodes.map(node => [node.id, node]));
      const linksCopy: LayoutLink[] = [];
      const linkedNodeIds = new Set<string>();
      let maxDist = 1;
      for (const link of props.links) {
        const source = nodeById.get(link.source);
        const target = nodeById.get(link.target);
        if (!source || !target || !Number.isFinite(link.snp_distance) || link.snp_distance < 0) continue;
        linksCopy.push({ ...link, source, target });
        linkedNodeIds.add(source.id);
        linkedNodeIds.add(target.id);
        maxDist = Math.max(maxDist, link.snp_distance);
      }

      const outerRadius = Math.max(240, 60 * Math.sqrt(linkedNodeIds.size) + 100);
      const activeSimulation = d3.forceSimulation<LayoutNode, LayoutLink>(layoutNodes)
        .stop()
        .force("link", d3.forceLink<LayoutNode, LayoutLink>(linksCopy).id(node => node.id).distance(80))
        .force("charge", d3.forceManyBody<LayoutNode>().strength(-110).distanceMax(200))
        // Centre linked samples and keep isolated samples outside the linked groups.
        // Avoid recentering the whole graph, which would pull the linked groups off centre.
        .force("radial", linkedNodeIds.size > 0
          ? d3.forceRadial<LayoutNode>(node => linkedNodeIds.has(node.id) ? 0 : outerRadius, 0, 0)
            .strength(node => linkedNodeIds.has(node.id) ? 0.04 : 0.12)
          : null)
        .force("center", linkedNodeIds.size === 0 ? d3.forceCenter<LayoutNode>(0, 0) : null)
        .force("collision", d3.forceCollide<LayoutNode>(40).iterations(2));
      simulation = activeSimulation;

      const g = svg.append("g");
      contentGroup = g;

      zoomBehavior = d3.zoom<SVGSVGElement, unknown>()
        .extent([[0, 0], [viewportWidth, viewportHeight]])
        .scaleExtent([0.2, 5])
        .on("start", (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
          if (!disposed && simulation === activeSimulation && event.sourceEvent) cameraMode = "manual";
        })
        .on("zoom", (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
          if (!disposed && simulation === activeSimulation) g.attr("transform", event.transform.toString());
        });
      svg.call(zoomBehavior);
      // Rebinding zoom otherwise retains the SVG's previous internal transform.
      applyTransform(d3.zoomIdentity);

      // Links
      const link = g.append("g")
        .selectAll<SVGLineElement, LayoutLink>("line")
        .data(linksCopy)
        .join("line")
        .attr("stroke", "#94a3b8")
        .attr("stroke-width", 1.5)
        .attr("stroke-opacity", d => 0.3 + 0.7 * (1 - d.snp_distance / maxDist));

      // Link labels
      const linkLabel = g.append("g")
        .selectAll<SVGTextElement, LayoutLink>("text")
        .data(linksCopy)
        .join("text")
        .attr("font-size", 9)
        .attr("fill", "#64748b")
        .attr("text-anchor", "middle")
        .text(d => d.snp_distance);

      // Node circles
      const node = g.append("g")
        .selectAll<SVGCircleElement, LayoutNode>("circle")
        .data(layoutNodes)
        .join("circle")
        .attr("r", 12)
        .attr("fill", d => clusterColor(d.cluster))
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5)
        .call(d3.drag<SVGCircleElement, LayoutNode>()
            .on("start", (event, d) => {
              if (disposed || simulation !== activeSimulation) return;
              cameraMode = "manual";
              if (!event.active) activeSimulation.alphaTarget(0.3).restart();
              d.fx = d.x; d.fy = d.y;
            })
            .on("drag", (event, d) => {
              if (disposed || simulation !== activeSimulation) return;
              d.fx = event.x; d.fy = event.y;
            })
            .on("end", (event, d) => {
              if (disposed || simulation !== activeSimulation) return;
              if (!event.active) activeSimulation.alphaTarget(0);
              d.fx = null; d.fy = null;
            })
        );

      // Node labels
      const label = g.append("g")
        .selectAll<SVGTextElement, LayoutNode>("text")
        .data(layoutNodes)
        .join("text")
        .attr("font-size", 10)
        .attr("fill", "#1e293b")
        .attr("dy", "0.35em")
        .attr("x", 11)
        .text(d => d.id.length > 20 ? d.id.slice(0, 18) + "…" : d.id);

      // Tooltip
      node.append("title").text(d => `${d.id} (Cluster ${d.cluster})`);

      function renderPositions(): void {
        if (disposed || simulation !== activeSimulation) return;
        link
          .attr("x1", d => d.source.x!)
          .attr("y1", d => d.source.y!)
          .attr("x2", d => d.target.x!)
          .attr("y2", d => d.target.y!);
        linkLabel
          .attr("x", d => (d.source.x! + d.target.x!) / 2)
          .attr("y", d => (d.source.y! + d.target.y!) / 2);
        node.attr("cx", d => d.x!).attr("cy", d => d.y!);
        label.attr("x", d => d.x! + 11).attr("y", d => d.y!);
        canFit.value = viewportVisible && hasFinitePositions();
      }

      activeSimulation
        .on("tick", renderPositions)
        .on("end", () => {
          if (!disposed && simulation === activeSimulation && cameraMode === "fit") fitGraph();
        });
      renderPositions();
      fitGraph();
      activeSimulation.restart();
    }

    function updateViewport(): void {
      if (disposed || !viewportRef.value || !svgRef.value) return;
      const width = viewportRef.value.clientWidth;
      const height = viewportRef.value.clientHeight;
      const wasVisible = viewportVisible;
      viewportVisible = width > 0 && height > 0;
      if (!viewportVisible) {
        canFit.value = false;
        simulation?.stop();
        return;
      }

      const previousWidth = viewportWidth;
      const previousHeight = viewportHeight;
      const transform = d3.zoomTransform(svgRef.value);
      viewportWidth = width;
      viewportHeight = height;
      if (needsBuild) {
        buildGraph();
        return;
      }

      canFit.value = hasFinitePositions();
      zoomBehavior?.extent([[0, 0], [width, height]]);
      if (!wasVisible || width !== previousWidth || height !== previousHeight) {
        if (cameraMode === "fit") {
          fitGraph();
        } else {
          const [centreX, centreY] = transform.invert([previousWidth / 2, previousHeight / 2]);
          applyTransform(d3.zoomIdentity.translate(
            width / 2 - transform.k * centreX,
            height / 2 - transform.k * centreY,
          ).scale(transform.k));
        }
      }
      if (!wasVisible && simulation && simulation.alpha() >= simulation.alphaMin()) simulation.restart();
    }

    function scheduleResize(): void {
      if (disposed || resizeFrame !== null) return;
      resizeFrame = requestAnimationFrame(() => {
        resizeFrame = null;
        updateViewport();
      });
    }

    onMounted(() => {
      resizeObserver = new ResizeObserver(scheduleResize);
      if (viewportRef.value) resizeObserver.observe(viewportRef.value);
      updateViewport();
    });

    onBeforeUnmount(() => {
      disposed = true;
      resizeObserver?.disconnect();
      resizeObserver = null;
      if (resizeFrame !== null) cancelAnimationFrame(resizeFrame);
      resizeFrame = null;
      disposeGraph();
    });

    watch([() => props.nodes, () => props.links], () => {
      if (disposed) return;
      disposeGraph();
      needsBuild = true;
      cameraMode = "fit";
      updateViewport();
    }, { flush: "post" });

    return { viewportRef, svgRef, canFit, fitGraph, clusterIds, clusterColor };
  },
});
</script>
