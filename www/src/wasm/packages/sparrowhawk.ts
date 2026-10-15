import { loadWebWasmPackage, WasmPackage } from "../loader";

// The assembler: wasm64 when the browser can run it, else the wasm-pack wasm32 build.

type SparrowhawkModule = typeof import("@/pkg/index.js");

export const sparrowhawkPackage: WasmPackage<SparrowhawkModule> = {
    id: "sparrowhawk",
    label: "Sparrowhawk",
    loadWasm64: () => loadWebWasmPackage<SparrowhawkModule>(
        new URL("/pkg_wasm64/sparrowhawk/index.js", globalThis.location.origin),
    ),
    async loadWasm32() {
        const [api, exports] = await Promise.all([
            import(/* webpackChunkName: "sparrowhawk-wasm32-fallback" */ "@/pkg/index.js"),
            import("@/pkg/index_bg.wasm"),
        ]);
        return { api, memory: exports.memory };
    },
};
