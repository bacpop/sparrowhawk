export {};

declare global {
    interface ElectronSaveRequest {
        suggestedName: string;
        data: ArrayBuffer;
        mimeType?: string;
    }

    interface ElectronSaveResponse {
        canceled: boolean;
        filePath?: string;
    }

    interface Window {
        electronAPI?: {
            saveFile(request: ElectronSaveRequest): Promise<ElectronSaveResponse>;
            readTaxonomicIndex(): Promise<ArrayBuffer>;
        };
    }
}
