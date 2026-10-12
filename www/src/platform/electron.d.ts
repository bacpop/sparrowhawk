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

    interface ElectronStreamSaveStartResponse {
        canceled: boolean;
        streamId?: string;
    }

    interface Window {
        electronAPI?: {
            saveFile(request: ElectronSaveRequest): Promise<ElectronSaveResponse>;
            beginAlignmentSave(suggestedName: string): Promise<ElectronStreamSaveStartResponse>;
            writeAlignmentSave(streamId: string, data: ArrayBuffer): Promise<void>;
            finishAlignmentSave(streamId: string): Promise<void>;
            abortAlignmentSave(streamId: string): Promise<void>;
            readTaxonomicIndex(): Promise<ArrayBuffer>;
        };
    }
}
