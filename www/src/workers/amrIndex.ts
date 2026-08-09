import { gunzipSync } from "fflate";

export const AMR_INDEX_FILE_NAME = "amrfinderplus_2026-05-15.1_dna_k23_amr-stress-virulence.amridx.gz";
export const AMR_INDEX_URL = `/${AMR_INDEX_FILE_NAME}`;

// eslint-disable-next-line @typescript-eslint/no-explicit-any
export async function createAmrDetector(wasm: any): Promise<any> {
    const response = await fetch(AMR_INDEX_URL);
    if (!response.ok) {
        throw new Error("index");
    }
    const compressed = new Uint8Array(await response.arrayBuffer());
    const indexBytes = gunzipSync(compressed);
    return new wasm.AmrDetector(indexBytes);
}
