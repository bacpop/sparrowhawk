import { loadWebWasmPackage, WasmPackage } from "../loader";


// Here we define the instance of WasmPackage per each of the modules we want to try to load as wasm64

type SkaModule = typeof import("@/pkg_ska/index.js");

export const skaPackage: WasmPackage<SkaModule> = {
    id: "ska",
    label: "SKA",
    loadWasm64: () => loadWebWasmPackage<SkaModule>(
        new URL("/pkg_wasm64/ska/index.js", globalThis.location.origin),
    ),
    async loadWasm32() {
        // Literal imports keep the existing automatically compiled fallback.
        const [api, exports] = await Promise.all([
            import(/* webpackChunkName: "ska-wasm32-fallback" */ "@/pkg_ska/index.js"),
            import("@/pkg_ska/index_bg.wasm"),
        ]);
        return { api, memory: exports.memory };
    },
};
