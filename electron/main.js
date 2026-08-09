const {app, BrowserWindow, dialog, ipcMain, protocol} = require("electron");
const fs = require("fs/promises");
const path = require("path");

const TAXONOMIC_INDEX_ASSET = "inverted_k_17_ss_50.ski";
const APP_SCHEME = "app";
const RENDERER_ROOT = path.join(__dirname, "renderer");

// Serve the packaged renderer from app:// instead of file://, so fetch() and workers
// resolve against a real, non-null origin.
protocol.registerSchemesAsPrivileged([
    {scheme: APP_SCHEME, privileges: {standard: true, secure: true, supportFetchAPI: true, stream: true}},
]);

// Chromium gates WebGPU behind Vulkan on GNU/Linux; without this the embedder falls back
// to CPU. Needs Electron >= 43 (older Chromium also required --enable-unsafe-webgpu).
if (process.platform === "linux") {
    app.commandLine.appendSwitch("enable-features", "Vulkan,DefaultANGLEVulkan,VulkanFromANGLE");
}

// application/wasm is required for WebAssembly.instantiateStreaming.
const MIME_TYPES = {
    ".html": "text/html",
    ".js": "text/javascript",
    ".css": "text/css",
    ".json": "application/json",
    ".wasm": "application/wasm",
    ".png": "image/png",
    ".svg": "image/svg+xml",
    ".ico": "image/x-icon",
    ".woff2": "font/woff2",
    ".webmanifest": "application/manifest+json",
    ".map": "application/json",
};

function registerAppProtocol() {
    // fs.readFile reads transparently from inside app.asar.
    protocol.handle(APP_SCHEME, async (request) => {
        const pathname = decodeURIComponent(new URL(request.url).pathname);
        const relative = pathname.replace(/^\/+/, "") || "index.html";
        const filePath = path.normalize(path.join(RENDERER_ROOT, relative));
        if (filePath !== RENDERER_ROOT && !filePath.startsWith(RENDERER_ROOT + path.sep)) {
            return new Response("Forbidden", {status: 403});
        }
        try {
            const data = await fs.readFile(filePath);
            const type = MIME_TYPES[path.extname(filePath).toLowerCase()] ?? "application/octet-stream";
            return new Response(data, {headers: {"content-type": type}});
        } catch {
            return new Response("Not found", {status: 404});
        }
    });
}

async function createMainWindow() {
    const win = new BrowserWindow({
        width: 1440,
        height: 960,
        minWidth: 1100,
        minHeight: 760,
        autoHideMenuBar: true,
        webPreferences: {
            contextIsolation: true,
            nodeIntegration: false,
            preload: path.join(__dirname, "preload.js"),
        },
    });

    if (process.env.ELECTRON_RENDERER_URL) {
        await win.loadURL(process.env.ELECTRON_RENDERER_URL);
        win.webContents.openDevTools({mode: "detach"});
        return;
    }

    await win.loadURL(`${APP_SCHEME}://bundle/index.html`);
}

ipcMain.handle("dialog:save-file", async (_event, {suggestedName, data}) => {
    const result = await dialog.showSaveDialog({
        defaultPath: suggestedName,
    });

    if (result.canceled || !result.filePath) {
        return {canceled: true};
    }

    await fs.writeFile(result.filePath, Buffer.from(new Uint8Array(data)));
    return {canceled: false, filePath: result.filePath};
});

ipcMain.handle("asset:read-taxonomic-index", async () => {
    const assetPath = process.env.ELECTRON_RENDERER_URL
        ? path.resolve(__dirname, "..", "www", "public", TAXONOMIC_INDEX_ASSET)
        : path.join(RENDERER_ROOT, TAXONOMIC_INDEX_ASSET);
    const bytes = await fs.readFile(assetPath);
    return Uint8Array.from(bytes).buffer;
});

app.whenReady().then(async () => {
    registerAppProtocol();
    await createMainWindow();

    app.on("activate", async () => {
        if (BrowserWindow.getAllWindows().length === 0) {
            await createMainWindow();
        }
    });
});

app.on("window-all-closed", () => {
    if (process.platform !== "darwin") {
        app.quit();
    }
});
