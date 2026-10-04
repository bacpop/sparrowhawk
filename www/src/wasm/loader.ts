import type { WasmRuntimeNotice } from "../types";


// This code in this folder allows to catch if a module compiled for wasm64 cannot be loaded, and tries to load a wasm32 one.
export interface WasmInstance<T> {
    api: T;
    memory?: WebAssembly.Memory;
}

export interface WasmPackage<T> {
    id: string;
    label: string;
    loadWasm32: () => Promise<WasmInstance<T>>;
    loadWasm64?: () => Promise<WasmInstance<T>>;
}

export interface WasmLoadResult<T> extends WasmInstance<T> {
    target: "wasm32" | "wasm64";
    fallbackReason?: string;
}


// Each worker has its own cache and WebAssembly instance.
const loads = new Map<string, Promise<WasmLoadResult<unknown>>>();

function errorMessage(error: unknown): string {
    return error instanceof Error ? error.message : String(error);
}

async function initialise<T>(spec: WasmPackage<T>): Promise<WasmLoadResult<T>> {
    let fallbackReason: string | undefined;
    if (spec.loadWasm64) {
        try {
            return { ...await spec.loadWasm64(), target: "wasm64" };
        } catch (error) {
            fallbackReason = errorMessage(error);
            console.warn(`[wasm] ${spec.label}: wasm64 loading failed; trying wasm32`, error);
        }
    }
    try {
        return { ...await spec.loadWasm32(), target: "wasm32", fallbackReason };
    } catch (error) {
        throw new Error(`${spec.label} could not be loaded: ${fallbackReason === undefined
            ? "" : `wasm64: ${fallbackReason}; `}wasm32: ${errorMessage(error)}`);
    }
}


export function loadWasmPackage<T>(
    spec: WasmPackage<T>,
    report?: (notice: WasmRuntimeNotice) => void,
): Promise<WasmLoadResult<T>> {
    let pending = loads.get(spec.id);
    if (!pending) {
        pending = initialise(spec);
        loads.set(spec.id, pending);
        // Allow a subsequent request to retry after both variants failed.
        void pending.catch(() => {
            if (loads.get(spec.id) === pending) loads.delete(spec.id);
        });
    }
    return (pending as Promise<WasmLoadResult<T>>).then((result) => {
        report?.({
            type: "wasm-runtime",
            moduleId: spec.id,
            label: spec.label,
            target: result.target,
            fallbackReason: result.fallbackReason,
        });
        return result;
    });
}


/** Load wasm-bindgen's `web` output directly, outside webpack's WASM pipeline. */
export async function loadWebWasmPackage<T>(moduleUrl: string | URL): Promise<WasmInstance<T>> {
    const url = new URL(moduleUrl, globalThis.location.href);
    const api = await import(/* webpackIgnore: true */ url.href);
    if (typeof api.default !== "function") {
        throw new Error(`Missing wasm-bindgen initializer in ${url.href}`);
    }
    const exports = await api.default({ module_or_path: new URL("index_bg.wasm", url) });
    return { api: api as T, memory: exports.memory };
}
