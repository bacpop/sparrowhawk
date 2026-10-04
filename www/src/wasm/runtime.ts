import type { WasmRuntimeNotice } from "../types";


// This defines a listener of events to see if fails or not

export function isWasmRuntimeNotice(value: unknown): value is WasmRuntimeNotice {
    if (typeof value !== "object" || value === null) return false;
    const notice = value as Partial<WasmRuntimeNotice>;
    return notice.type === "wasm-runtime" && typeof notice.moduleId === "string"
        && typeof notice.label === "string"
        && (notice.target === "wasm32" || notice.target === "wasm64")
        && (notice.fallbackReason === undefined || typeof notice.fallbackReason === "string");
}

/** Install immediately after worker creation, before job handlers are assigned. */
export function observeWasmRuntime(
    worker: Worker,
    record: (notice: WasmRuntimeNotice) => void,
): () => void {
    const listener = (event: MessageEvent<unknown>) => {
        if (!isWasmRuntimeNotice(event.data)) return;
        // Job handlers otherwise mistake startup notices for job results.
        event.stopImmediatePropagation();
        record(event.data);
    };
    worker.addEventListener("message", listener);
    return () => worker.removeEventListener("message", listener);
}

export function retainWasmRuntimeStatus(
    previous: WasmRuntimeNotice | undefined,
    notice: WasmRuntimeNotice,
): WasmRuntimeNotice {
    // Pools can contain workers using different variants. Keep evidence of any fallback.
    return previous?.fallbackReason !== undefined && notice.fallbackReason === undefined
        ? previous : { ...notice };
}
