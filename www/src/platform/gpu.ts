/**
 * WebGPU devices detection and probing
 */

import type { GpuAdapterInfo } from "@/types";

export type { GpuAdapterInfo };

export function isWebGpuAvailable(): boolean {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const gpu = (navigator as any)?.gpu;
    const secure = typeof globalThis.isSecureContext === "boolean"
        ? globalThis.isSecureContext
        : true;
    return Boolean(gpu) && secure;
}


// Via WebGL, as the info from JS might be short
export function getWebGLRendererLabel(): string {
    try {
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        let canvas: any;
        if (typeof OffscreenCanvas !== "undefined") {
            canvas = new OffscreenCanvas(1, 1);
        } else if (typeof document !== "undefined") {
            canvas = document.createElement("canvas");
        } else {
            return '';
        }

        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const gl: any = canvas.getContext('webgl2') || canvas.getContext('webgl');
        if (!gl) return '';
        
        // Firefox deprecated WEBGL_debug_renderer_info and now sanitises RENDERER to the same
        // string. Chrome still masks RENDERER, so we check both just in case.
        let renderer: string = gl.getParameter(gl.RENDERER) || '';
        if (!renderer || /^(webkit\s+)?webgl$|^mozilla$/i.test(renderer.trim())) {
            const ext = gl.getExtension('WEBGL_debug_renderer_info');
            if (ext) renderer = gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) || renderer;
        }
        if (!renderer) return '';

        // Extract the GPU model from the ANGLE wrapper when present:
        //   "ANGLE (Vendor, <model>, Backend)" → "<model>"
        const m = renderer.match(/^ANGLE\s*\([^,]+,\s*(.+),\s*[^,]+\)$/);
        if (m) renderer = m[1];

        // Strip backend-description suffixes that describe rendering API, not the GPU:
        renderer = renderer
            .replace(/^ANGLE\s+Metal\s+Renderer:\s*/i, '')
            .replace(/\s+Direct3D\d+(\s+vs_\d+_\d+\s+ps_\d+_\d+)?$/i, '')
            .trim();

        return renderer;
    } catch {
        return '';
    }
}

// JS API, uses the other methods as well to fill info
export async function listGpuAdapters(): Promise<GpuAdapterInfo[]> {
    if (!isWebGpuAvailable()) {
        console.warn("[Sparrowhawk] WebGPU not available (no navigator.gpu, or insecure context)");
        return [];
    }
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const gpu = (navigator as any).gpu;

    // Lazy: only adapters the browser refused to name need this, so Chrome never creates a
    // WebGL context or touches the deprecated extension.
    let cachedRenderer: string | undefined;
    const webglRenderer = () => (cachedRenderer ??= getWebGLRendererLabel());

    // Pass 1: collect (index, hwLabel) for every power-preference that returns an adapter.
    interface AdapterEntry { index: number; hwLabel: string; identified: boolean; fingerprint: string }
    const entries: AdapterEntry[] = [];

    const prefs: [number, string | undefined][] = [
        [0, undefined],
        [1, 'high-performance'],
        [2, 'low-power'],
    ];

    for (const [idx, powerPreference] of prefs) {
        try {
            // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const adapter: any = await gpu.requestAdapter(
                powerPreference !== undefined ? { powerPreference } : undefined
            );
            if (adapter) {
                const info = adapter.info;
                const description = (info?.description || '').trim();
                const vendorArch = (info?.vendor && info?.architecture)
                    ? `${info.vendor} (${info.architecture})`
                    : '';
                const vendorDevice = info?.vendor
                    ? `${info.vendor}${info?.device ? ` ${info.device}` : ''}`
                    : (info?.device || '');
                const hwLabel = description || vendorArch || vendorDevice || webglRenderer() || 'GPU';
                // The WebGL fallback is one global string, so it cannot tell adapters apart.
                const identified = Boolean(description || vendorArch || vendorDevice);
                // limits and features survive info redaction, so they are what distinguishes
                // adapters when the name cannot.
                const fingerprint = [
                    adapter.limits?.maxComputeWorkgroupStorageSize,
                    adapter.limits?.maxBufferSize,
                    [...(adapter.features ?? [])].sort().join(","),
                ].join("|");
                entries.push({ index: idx, hwLabel, identified, fingerprint });
            }
        } catch (e) {
            console.warn(`[Sparrowhawk] requestAdapter(${powerPreference ?? 'default'}) failed:`, e);
        }
    }

    // Only two outcomes are reachable: cubecl maps IntegratedGpu to low power and everything
    // else to high performance, so the "default" probe always duplicates high-performance.
    const hp = entries.find(e => e.index === 1) ?? entries.find(e => e.index === 0);
    const lp = entries.find(e => e.index === 2);

    // Identity comes from the fingerprint, not the label: two probes can return one device
    // under different names, or two devices under none.
    const found: { entry: AdapterEntry; index: number; role: string }[] = [];
    if (hp) found.push({ entry: hp, index: 1, role: "High performance" });
    if (lp && !(hp && hp.fingerprint === lp.fingerprint)) {
        found.push({ entry: lp, index: 2, role: "Low power" });
    }

    // The WebGL renderer is one global string, so it can only name an adapter when there is a
    // single one, if more, there is uncertainty
    const sole = found.length === 1 ? webglRenderer() : "";
    const results: GpuAdapterInfo[] = found.map(({ entry, index, role }) => ({
        index,
        name: entry.identified ? entry.hwLabel : sole || role,
        identified: entry.identified,
    }));

    console.log("[Sparrowhawk] Detected GPU adapters:", results);
    return results;
}
