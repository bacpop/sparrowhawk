declare module "@/pkg_mandrake/index" {
  export function clusterEmbedding(embedding: Float64Array): Int32Array;
  export function sketchKmerLengths(file: File): Uint32Array;

  export class MandrakeProgress {
    readonly completed: number;
    readonly maximum: number;
    readonly eq: number;
    readonly complete: boolean;
  }

  export class MandrakeDistanceProgress {
    readonly completed: number;
    readonly maximum: number;
    readonly complete: boolean;
  }

  export class MandrakeOperation {
    static fromAlignmentFile(
      file: File,
      mode: string,
      value: number,
      perplexity: number,
      maxUpdates: number,
      repulsionSamples: number,
      learningRate: number,
      initialExaggeration: boolean,
    ): MandrakeOperation;
    static fromAccessoryFile(
      file: File,
      mode: string,
      value: number,
      perplexity: number,
      maxUpdates: number,
      repulsionSamples: number,
      learningRate: number,
      initialExaggeration: boolean,
    ): MandrakeOperation;
    static fromSketchFiles(
      metadataFile: File,
      dataFile: File,
      mode: string,
      value: number,
      perplexity: number,
      maxUpdates: number,
      repulsionSamples: number,
      learningRate: number,
      initialExaggeration: boolean,
      distanceKind: string,
      jaccardKmer: number,
    ): MandrakeOperation;
    advanceDistances(rowBudget: number): MandrakeDistanceProgress;
    beginEmbedding(): void;
    advance(roundBudget: number): MandrakeProgress;
    embedding(): Float64Array;
    names(): string;
    sample_count(): number;
    is_complete(): boolean;
  }
}
