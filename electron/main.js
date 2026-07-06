const {app, BrowserWindow, dialog, ipcMain} = require("electron");
const fs = require("fs/promises");
const path = require("path");

const TAXONOMIC_INDEX_ASSET = "inverted_k_17_ss_50.ski";

function getRendererEntry() {
    if (process.env.ELECTRON_RENDERER_URL) {
        return process.env.ELECTRON_RENDERER_URL;
    }

    return path.join(__dirname, "renderer", "index.html");
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

    const entry = getRendererEntry();
    if (/^https?:\/\//.test(entry)) {
        await win.loadURL(entry);
        win.webContents.openDevTools({mode: "detach"});
        return;
    }

    await win.loadFile(entry);
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
        : path.join(__dirname, "renderer", TAXONOMIC_INDEX_ASSET);
    const bytes = await fs.readFile(assetPath);
    return Uint8Array.from(bytes).buffer;
});

app.whenReady().then(async () => {
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
