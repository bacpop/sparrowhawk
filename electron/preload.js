const {contextBridge, ipcRenderer} = require("electron");

contextBridge.exposeInMainWorld("electronAPI", {
    saveFile(request) {
        return ipcRenderer.invoke("dialog:save-file", request);
    },
    readTaxonomicIndex() {
        return ipcRenderer.invoke("asset:read-taxonomic-index");
    },
});
