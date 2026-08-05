import { gunzipSync } from "fflate";
import { loadAssetBytes } from "@/platform/files";

/** Staged by `rust/esm-bridge/scripts/stage_weights.sh`. f16, 13.7 MB gzipped. */
export const ESM_WEIGHTS_FILE_NAME = "esm2_t6_8M_UR50D.bpkz";
export const ESM_ENCODER_FILE_NAME = "esm2_umap_encoder.bpkz";

let cached: Uint8Array | null = null;
let cachedEncoder: Uint8Array | undefined;

export async function loadEsmWeights(): Promise<Uint8Array> {
    if (cached) return cached;
    cached = gunzipSync(await loadAssetBytes(ESM_WEIGHTS_FILE_NAME));
    return cached;
}

export async function loadEsmEncoder(): Promise<Uint8Array> {
    if (cachedEncoder !== undefined) return cachedEncoder;
    cachedEncoder = gunzipSync(await loadAssetBytes(ESM_ENCODER_FILE_NAME));
    return cachedEncoder;
}
