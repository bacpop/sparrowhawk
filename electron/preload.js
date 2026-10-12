const {contextBridge, ipcRenderer} = require("electron");

contextBridge.exposeInMainWorld("electronAPI", {
    saveFile(request) {
        return ipcRenderer.invoke("dialog:save-file", request);
    },
    beginAlignmentSave(suggestedName) {
        return ipcRenderer.invoke("dialog:begin-alignment-save", {suggestedName});
    },
    writeAlignmentSave(streamId, data) {
        return ipcRenderer.invoke("dialog:write-alignment-save", {streamId, data});
    },
    finishAlignmentSave(streamId) {
        return ipcRenderer.invoke("dialog:finish-alignment-save", {streamId});
    },
    abortAlignmentSave(streamId) {
        return ipcRenderer.invoke("dialog:abort-alignment-save", {streamId});
    },
    readTaxonomicIndex() {
        return ipcRenderer.invoke("asset:read-taxonomic-index");
    },
});
