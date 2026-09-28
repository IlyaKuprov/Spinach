# kernel/includes/parallel_profiler_start.m

- Signature: `(script file)`

## Purpose

Start parallel-pool byte counting and timing around a profiled parallel stage.

## Behaviour

On non-worker nodes, the include calls `ticBytes(gcp)`; it calls `tic()` on every node. When running on a non-worker node with `'dafuq'` enabled in `spin_system.sys.enable`, it creates a `PoolProfiler` object in `parProfiler`.

## Source documentation

https://spindynamics.org/wiki/index.php?title=parallel_profiler_start.m
