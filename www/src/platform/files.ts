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

    // Against the origin, not the href: a worker's href can be an opaque `blob:` URL: for the embeddings
    if (location.origin && location.origin !== "null") {
        return new URL("/" + name, location.origin).toString();
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

function describeCause(cause: unknown): string {
    return cause instanceof Error ? cause.message : String(cause);
}

/** Resolve, fetch, and check the status. Shared by both asset loaders below. */
async function fetchAsset(name: string): Promise<Response> {
    const url = getAssetUrl(name);

    // Proper logging with the reason why we cannot load something
    let response: Response;
    try {
        response = await fetch(url);
    } catch (cause) {
        throw new Error(`Could not fetch asset ${name} from ${url}: ${describeCause(cause)}`);
    }

    if (!response.ok) {
        throw new Error(
            `Failed to load asset ${name} from ${url}: ${response.status} ${response.statusText}`,
        );
    }

    return response;
}

// A helper to catch errors while reading a file that is being loaded
async function readAssetBody<T>(
    name: string,
    response: Response,
    read: (response: Response) => Promise<T>,
): Promise<T> {
    try {
        return await read(response);
    } catch (cause) {
        const expected = response.headers.get("content-length");
        throw new Error(
            `Asset ${name} failed while reading its body` +
            (expected ? ` (expected ${expected} bytes)` : "") +
            `: ${describeCause(cause)}`,
        );
    }
}

// The inverted index is loaded in a particular way on the Electron version to avoid issues
export async function loadAssetBlob(name: string): Promise<Blob> {
    if (
        typeof window !== "undefined" &&
        typeof window.electronAPI?.readTaxonomicIndex === "function" &&
        name === "inverted_k_17_ss_50.ski"
    ) {
        const data = await window.electronAPI.readTaxonomicIndex();
        return new Blob([data]);
    }

    const response = await fetchAsset(name);
    return readAssetBody(name, response, (r) => r.blob());
}

// This fixes loading large files in Chromium, directly as bytes instead of as a blob
export async function loadAssetBytes(name: string): Promise<Uint8Array> {
    const response = await fetchAsset(name);
    return new Uint8Array(await readAssetBody(name, response, (r) => r.arrayBuffer()));
}
