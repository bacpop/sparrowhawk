import * as d3 from "d3";

// Tableau10 without its grey (#bab0ac) and greyish brown (#9c755f), so clusters never look like singletons.
const BASE_CLUSTER_PALETTE: readonly string[] =
  d3.schemeTableau10.filter(color => color !== "#bab0ac" && color !== "#9c755f");

export const SINGLETON_COLOR = "#808080";
export const SINGLETON_LABEL = "Singleton";

export function transmissionClusterColor(clusterId: number): string {
  if (clusterId >= 1 && clusterId <= BASE_CLUSTER_PALETTE.length) return BASE_CLUSTER_PALETTE[clusterId - 1];
  // Golden-angle hue spacing: deterministic, unbounded, and always saturated (never grey).
  const hue = (((clusterId * 137.508) % 360) + 360) % 360;
  return d3.hsl(hue, 0.72, 0.52).formatHex();
}

type TaxoniumRgb = [number, number, number];

function hexToRgb(hex: string): TaxoniumRgb {
  return [1, 3, 5].map(i => Number.parseInt(hex.slice(i, i + 2), 16)) as TaxoniumRgb;
}

export function buildTaxoniumClusterColorMapping(clusterIds: number[]): Record<string, TaxoniumRgb> {
  const mapping: Record<string, TaxoniumRgb> = { [SINGLETON_LABEL]: hexToRgb(SINGLETON_COLOR) };
  for (const clusterId of new Set(clusterIds)) {
    mapping[`Cluster ${clusterId}`] = hexToRgb(transmissionClusterColor(clusterId));
  }
  return mapping;
}
