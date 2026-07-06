function hasElectronSave(): boolean {
    return typeof window !== "undefined" && typeof window.electronAPI?.saveFile === "function";
}

function toUint8Array(content: Uint8Array | ArrayBuffer): Uint8Array {
    return content instanceof Uint8Array ? content : new Uint8Array(content);
}

function toArrayBuffer(content: Uint8Array): ArrayBuffer {
    const out = new Uint8Array(content.byteLength);
    out.set(content);
    return out.buffer;
}

function triggerDownload(blob: Blob, filename: string): void {
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = filename;
    anchor.style.display = "none";
    document.body.appendChild(anchor);
    anchor.click();
    document.body.removeChild(anchor);
    URL.revokeObjectURL(url);
}

function getAssetUrl(name: string): string {
    const location = globalThis.location;
    const isWorkerContext = typeof window === "undefined";

    if (!location) {
        return "/" + name;
    }

    if (location.protocol === "file:") {
        const relativePath = (isWorkerContext ? "../" : "./") + name;
        return new URL(relativePath, location.href).toString();
    }

    return new URL("/" + name, location.href).toString();
}

export async function saveTextFile(content: string, filename: string, mimeType = "text/plain"): Promise<void> {
    if (hasElectronSave()) {
        const data = new TextEncoder().encode(content);
        await window.electronAPI!.saveFile({
            suggestedName: filename,
            data: toArrayBuffer(data),
            mimeType,
        });
        return;
    }

    triggerDownload(new Blob([content], {type: `${mimeType};charset=utf-8`}), filename);
}

export async function saveBinaryFile(
    content: Uint8Array | ArrayBuffer,
    filename: string,
    mimeType = "application/octet-stream"
): Promise<void> {
    const data = toUint8Array(content);

    if (hasElectronSave()) {
        await window.electronAPI!.saveFile({
            suggestedName: filename,
            data: toArrayBuffer(data),
            mimeType,
        });
        return;
    }

    triggerDownload(new Blob([toArrayBuffer(data)], {type: mimeType}), filename);
}

export async function loadAssetBlob(name: string): Promise<Blob> {
    if (
        typeof window !== "undefined" &&
        typeof window.electronAPI?.readTaxonomicIndex === "function" &&
        name === "inverted_k_17_ss_50.ski"
    ) {
        const data = await window.electronAPI.readTaxonomicIndex();
        return new Blob([data]);
    }

    const response = await fetch(getAssetUrl(name));
    if (!response.ok) {
        throw new Error(`Failed to load asset ${name}: ${response.status}`);
    }

    return response.blob();
}
