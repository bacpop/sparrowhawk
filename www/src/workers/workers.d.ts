// Type declarations for worker modules processed by worker-loader
declare module '@/workers/Assembler.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Mapper.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Sketcher.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Caller.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Depleter.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/AmrDetector.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Embedder.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}

declare module '@/workers/Mandrake.worker' {
    const WorkerFactory: new () => Worker;
    export default WorkerFactory;
}
