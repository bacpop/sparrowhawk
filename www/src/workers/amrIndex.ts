import { gunzipSync } from "fflate";
import type { AmrDetector } from "@/pkg_amr";

export const AMR_INDEX_FILE_NAME = "amrfinderplus_2026-05-15.1_dna_k23_amr-stress-virulence.amridx.gz";
export const AMR_INDEX_URL = `/${AMR_INDEX_FILE_NAME}`;

export async function createAmrDetector(wasm: typeof import("@/pkg_amr")): Promise<AmrDetector> {
    const response = await fetch(AMR_INDEX_URL);
    if (!response.ok) {
        throw new Error("index");
    }
    const compressed = new Uint8Array(await response.arrayBuffer());
    const indexBytes = gunzipSync(compressed);
    return new wasm.AmrDetector(indexBytes);
}
