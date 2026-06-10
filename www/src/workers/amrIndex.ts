import { gunzipSync } from "fflate";

export const AMR_INDEX_FILE_NAME = "amrfinderplus_2026-03-24.1_dna_k31_amr-stress-virulence.amridx.gz";
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
