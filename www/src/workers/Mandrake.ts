import WorkerMandrake from "@/workers/Mandrake.worker";

export interface MandrakeSettings {
  mode: "knn" | "threshold";
  value: number;
  perplexity: number;
  maxUpdates: number;
  repulsionSamples: number;
  learningRate: number;
  initialExaggeration: boolean;
  hdbscan: boolean;
  sketchDistance: "core" | "jaccard";
  jaccardKmer: number;
}

export interface MandrakeSketchFiles {
  metadata: File;
  data: File;
}

export interface MandrakeProgress {
  completed: number;
  maximum: number;
  eq: number;
  complete: boolean;
}

export interface MandrakeDistanceProgress {
  completed: number;
  maximum: number;
  complete: boolean;
}

export interface MandrakeResult {
  embedding: Float64Array;
  names: string[];
  completed: number;
  maximum: number;
  eq: number;
  hdbscanLabels: Int32Array | null;
  hdbscanError: string | null;
}

export type MandrakeUpdate =
  | { phase: "sketch" }
  | { phase: "distance"; progress: MandrakeDistanceProgress; names?: string[] }
  | { phase: "embedding"; progress: MandrakeProgress }
  | { phase: "clustering" }
  | { phase: "frame"; embedding: Float64Array; completed: number; maximum: number };

type UpdateHandler = (update: MandrakeUpdate) => void;

interface DistanceProgressMessage extends MandrakeDistanceProgress {
  type: "distance-progress";
  names: string;
}

interface SketchLoadingMessage {
  type: "sketch-loading";
}

interface SketchMetadataMessage {
  type: "sketch-metadata";
  kmers: number[];
}

interface EmbeddingProgressMessage extends MandrakeProgress {
  type: "embedding-progress";
}

interface FrameMessage {
  type: "frame";
  embedding: Float64Array;
  completed: number;
  maximum: number;
}

interface CompleteMessage {
  type: "complete";
  embedding: Float64Array;
  names: string;
  completed: number;
  maximum: number;
  eq: number;
  clusters?: Int32Array;
  hdbscanError?: string;
}

interface ClusteringMessage {
  type: "clustering";
}

interface ErrorMessage {
  type: "error";
  message: string;
}

type WorkerMessage =
  | SketchLoadingMessage
  | SketchMetadataMessage
  | DistanceProgressMessage
  | EmbeddingProgressMessage
  | ClusteringMessage
  | FrameMessage
  | CompleteMessage
  | ErrorMessage;

export class MandrakeRunner {
  private worker: Worker | null = null;
  private updateHandler: UpdateHandler | null = null;
  private resolve: ((result: MandrakeResult) => void) | null = null;
  private reject: ((error: Error) => void) | null = null;
  private distanceRowBudget = 1;
  private embeddingRoundBudget = 64;

  run(
    source: "alignment" | "accessory",
    file: File,
    settings: MandrakeSettings,
    onUpdate: UpdateHandler,
  ): Promise<MandrakeResult> {
    this.cancel();
    const worker = new WorkerMandrake();
    this.worker = worker;
    this.updateHandler = onUpdate;

    const result = new Promise<MandrakeResult>((resolve, reject) => {
      this.resolve = resolve;
      this.reject = reject;
    });
    worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
      this.handleMessage(event.data);
    };
    worker.onerror = (event: ErrorEvent) => {
      this.fail(new Error(event.message || "Mandrake worker failed"));
    };

    worker.postMessage({ type: "start", source, file, settings });
    return result;
  }

  runSketch(
    files: MandrakeSketchFiles,
    settings: MandrakeSettings,
    onUpdate: UpdateHandler,
  ): Promise<MandrakeResult> {
    this.cancel();
    const worker = new WorkerMandrake();
    this.worker = worker;
    this.updateHandler = onUpdate;

    const result = new Promise<MandrakeResult>((resolve, reject) => {
      this.resolve = resolve;
      this.reject = reject;
    });
    worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
      this.handleMessage(event.data);
    };
    worker.onerror = (event: ErrorEvent) => {
      this.fail(new Error(event.message || "Mandrake worker failed"));
    };

    worker.postMessage({ type: "start-sketch", metadata: files.metadata, data: files.data, settings });
    return result;
  }

  inspectSketchKmers(file: File): Promise<number[]> {
    const worker = new WorkerMandrake();
    return new Promise<number[]>((resolve, reject) => {
      worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
        if (event.data.type === "sketch-metadata") {
          worker.terminate();
          resolve(event.data.kmers);
        } else if (event.data.type === "error") {
          worker.terminate();
          reject(new Error(event.data.message));
        }
      };
      worker.onerror = (event: ErrorEvent) => {
        worker.terminate();
        reject(new Error(event.message || "Mandrake worker failed"));
      };
      worker.postMessage({ type: "inspect-sketch", file });
    });
  }

  cancel(): void {
    if (this.worker) {
      this.worker.terminate();
      this.worker = null;
    }
    if (this.reject) {
      this.reject(new Error("Mandrake operation cancelled"));
    }
    this.resolve = null;
    this.reject = null;
    this.updateHandler = null;
  }

  private handleMessage(message: WorkerMessage): void {
    if (message.type === "error") {
      this.fail(new Error(message.message));
      return;
    }

    if (message.type === "sketch-loading") {
      this.updateHandler?.({ phase: "sketch" });
      return;
    }

    if (message.type === "sketch-metadata") return;

    if (message.type === "distance-progress") {
      const progress: MandrakeDistanceProgress = {
        completed: message.completed,
        maximum: message.maximum,
        complete: message.complete,
      };
      this.updateHandler?.({
        phase: "distance",
        progress,
        names: message.names ? message.names.split("\n") : [],
      });
      this.distanceRowBudget = Math.max(1, Math.ceil(message.maximum / 100));
      if (message.complete) {
        this.worker?.postMessage({ type: "begin-embedding" });
      } else {
        this.worker?.postMessage({ type: "advance-distance", rowBudget: this.distanceRowBudget });
      }
      return;
    }

    if (message.type === "embedding-progress") {
      const progress: MandrakeProgress = {
        completed: message.completed,
        maximum: message.maximum,
        eq: message.eq,
        complete: message.complete,
      };
      this.updateHandler?.({ phase: "embedding", progress });
      if (!message.complete) {
        this.embeddingRoundBudget = Math.max(64, Math.ceil(message.maximum / 100));
        this.worker?.postMessage({ type: "advance-embedding", roundBudget: this.embeddingRoundBudget });
      }
      return;
    }

    if (message.type === "frame") {
      this.updateHandler?.({
        phase: "frame",
        embedding: message.embedding,
        completed: message.completed,
        maximum: message.maximum,
      });
      return;
    }

    if (message.type === "clustering") {
      this.updateHandler?.({ phase: "clustering" });
      return;
    }

    this.updateHandler?.({
      phase: "embedding",
      progress: {
        completed: message.completed,
        maximum: message.maximum,
        eq: message.eq,
        complete: true,
      },
    });
    const result: MandrakeResult = {
      embedding: message.embedding,
      names: message.names ? message.names.split("\n") : [],
      completed: message.completed,
      maximum: message.maximum,
      eq: message.eq,
      hdbscanLabels: message.clusters ?? null,
      hdbscanError: message.hdbscanError ?? null,
    };
    this.resolve?.(result);
    this.cleanupWorker();
  }

  private fail(error: Error): void {
    this.reject?.(error);
    this.cleanupWorker();
  }

  private cleanupWorker(): void {
    this.worker?.terminate();
    this.worker = null;
    this.resolve = null;
    this.reject = null;
    this.updateHandler = null;
  }
}
