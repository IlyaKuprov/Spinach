# kernel/includes/parallel_profiler_report.m

- Signature: `(script file)`

## Purpose

An include that writes a report of profiling around parallel stages. Invoke it just after a `parfor` or `spmd` stage for which `parallel_profiler_start` was called.

## Behaviour

On non-worker nodes, it reports the average worker-process data received and sent in MB and the elapsed parallel-stage time in seconds. If `'dafuq'` is enabled in `spin_system.sys.enable`, it drains the profiler log, saves `parpool_history` to a timestamped MAT file in `spin_system.sys.scratch`, and reports the filename.

## Source comment attribution

The source includes a poem attributed to Philip Larkin.

## Source documentation

https://spindynamics.org/wiki/index.php?title=parallel_profiler_report.m
