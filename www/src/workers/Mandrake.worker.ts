import type { MandrakeOperation } from "@/pkg_mandrake/index";

type MandrakeWasm = typeof import("@/pkg_mandrake/index");

interface Settings {
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

interface StartMessage {
  type: "start";
  source: "alignment" | "accessory";
  file: File;
  settings: Settings;
}

interface StartSketchMessage {
  type: "start-sketch";
  metadata: File;
  data: File;
  settings: Settings;
}

interface InspectSketchMessage {
  type: "inspect-sketch";
  file: File;
}

interface AdvanceDistanceMessage {
  type: "advance-distance";
  rowBudget: number;
}

interface BeginEmbeddingMessage {
  type: "begin-embedding";
}

interface AdvanceEmbeddingMessage {
  type: "advance-embedding";
  roundBudget: number;
}

interface ResetMessage {
  type: "reset";
}

type WorkerMessage =
  | StartMessage
  | StartSketchMessage
  | InspectSketchMessage
  | AdvanceDistanceMessage
  | BeginEmbeddingMessage
  | AdvanceEmbeddingMessage
  | ResetMessage;

let operation: MandrakeOperation | null = null;
let queue = Promise.resolve();
let wasmPromise: Promise<MandrakeWasm> | null = null;
let nextFrameUpdate = 0;
let frameInterval = 1;
let runHdbscan = false;
const MAX_FRAME_INTERVAL = 20_000;

function loadWasm(): Promise<MandrakeWasm> {
  if (!wasmPromise) wasmPromise = import("@/pkg_mandrake/index");
  return wasmPromise;
}

function reportError(error: unknown): void {
  const message = error instanceof Error ? error.message : String(error);
  self.postMessage({ type: "error", message });
}

function postDistanceProgress(): void {
  if (!operation) throw new Error("no Mandrake operation is active");
  const progress = operation.advanceDistances(0);
  self.postMessage({
    type: "distance-progress",
    completed: progress.completed,
    maximum: progress.maximum,
    complete: progress.complete,
    names: operation.names(),
  });
}

function postFrame(completed: number, maximum: number): void {
  if (!operation) throw new Error("no Mandrake operation is active");
  const embedding = operation.embedding();
  self.postMessage(
    { type: "frame", embedding, completed, maximum },
    { transfer: [embedding.buffer] },
  );
}

async function handle(message: WorkerMessage): Promise<void> {
  if (message.type === "reset") {
    operation = null;
    runHdbscan = false;
    self.postMessage({ type: "reset" });
    return;
  }

  if (message.type === "start") {
    const { MandrakeOperation } = await loadWasm();
    const settings = message.settings;
    runHdbscan = settings.hdbscan;
    operation = message.source === "alignment"
      ? MandrakeOperation.fromAlignmentFile(
          message.file,
          settings.mode,
          settings.value,
          settings.perplexity,
          settings.maxUpdates,
          settings.repulsionSamples,
          settings.learningRate,
          settings.initialExaggeration,
        )
      : MandrakeOperation.fromAccessoryFile(
          message.file,
          settings.mode,
          settings.value,
          settings.perplexity,
          settings.maxUpdates,
          settings.repulsionSamples,
          settings.learningRate,
          settings.initialExaggeration,
        );
    postDistanceProgress();
    return;
  }

  if (message.type === "inspect-sketch") {
    const wasm = await loadWasm();
    const kmers = wasm.sketchKmerLengths(message.file);
    self.postMessage({ type: "sketch-metadata", kmers });
    return;
  }

  if (message.type === "start-sketch") {
    const { MandrakeOperation } = await loadWasm();
    const settings = message.settings;
    runHdbscan = settings.hdbscan;
    self.postMessage({ type: "sketch-loading" });
    operation = MandrakeOperation.fromSketchFiles(
      message.metadata,
      message.data,
      settings.mode,
      settings.value,
      settings.perplexity,
      settings.maxUpdates,
      settings.repulsionSamples,
      settings.learningRate,
      settings.initialExaggeration,
      settings.sketchDistance,
      settings.jaccardKmer,
    );
    postDistanceProgress();
    return;
  }

  if (!operation) throw new Error("no Mandrake operation is active");

  if (message.type === "advance-distance") {
    const progress = operation.advanceDistances(message.rowBudget);
    self.postMessage({
      type: "distance-progress",
      completed: progress.completed,
      maximum: progress.maximum,
      complete: progress.complete,
      names: operation.names(),
    });
    return;
  }

  if (message.type === "begin-embedding") {
    operation.beginEmbedding();
    const progress = operation.advance(0);
    frameInterval = Math.max(1, Math.min(MAX_FRAME_INTERVAL, Math.ceil(progress.maximum / 20)));
    nextFrameUpdate = frameInterval;
    self.postMessage({
      type: "embedding-progress",
      completed: progress.completed,
      maximum: progress.maximum,
      eq: progress.eq,
      complete: progress.complete,
    });
    postFrame(progress.completed, progress.maximum);
    return;
  }

  if (message.type === "advance-embedding") {
    const progress = operation.advance(message.roundBudget);
    if (progress.complete) {
      const embedding = operation.embedding();
      let clusters: Int32Array | undefined;
      let hdbscanError: string | undefined;
      if (runHdbscan) {
        self.postMessage({ type: "clustering" });
        try {
          const wasm = await loadWasm();
          clusters = wasm.clusterEmbedding(embedding);
        } catch (error) {
          hdbscanError = error instanceof Error ? error.message : String(error);
        }
      }
      self.postMessage({
        type: "embedding-progress",
        completed: progress.completed,
        maximum: progress.maximum,
        eq: progress.eq,
        complete: true,
      });
      self.postMessage({
        type: "complete",
        completed: progress.completed,
        maximum: progress.maximum,
        eq: progress.eq,
        embedding,
        names: operation.names(),
        clusters,
        hdbscanError,
      }, {
        transfer: [
          embedding.buffer,
          ...(clusters ? [clusters.buffer] : []),
        ],
      });
      return;
    }

    self.postMessage({
      type: "embedding-progress",
      completed: progress.completed,
      maximum: progress.maximum,
      eq: progress.eq,
      complete: false,
    });
    if (progress.completed >= nextFrameUpdate) {
      while (nextFrameUpdate <= progress.completed) nextFrameUpdate += frameInterval;
      postFrame(progress.completed, progress.maximum);
    }
  }
}

self.onmessage = (event: MessageEvent<WorkerMessage>) => {
  queue = queue.then(() => handle(event.data)).catch(reportError);
};
