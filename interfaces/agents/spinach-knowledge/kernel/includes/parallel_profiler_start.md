# kernel/includes/parallel_profiler_start.m

- Signature: script include; no function output.

## Purpose

Start timing and, when requested, parallel-pool byte-counting/profiling immediately before a parallel stage. The source comment names `parfor` or `spmd` as the intended placement.

## Execution and data mapping

- On a non-worker node, call `ticBytes(gcp)` to start byte counting for the current parallel pool.
- Call `tic()` unconditionally on every node that executes the include, starting a local elapsed-time timer.
- Only on a non-worker node, and only when `'dafuq'` is a member of `spin_system.sys.enable`, construct `parallel.internal.profiling.PoolProfiler()` and assign it to `parProfiler`.
- This include starts measurements; it does not itself run the following parallel loop, stop the timers, report results, or modify the spin system or its operators.

## Source comments

The source also carries the unrelated Gauss-summation example: adding the integers 1 through 100 gives 5050. It is illustrative commentary, not profiler logic.

## References

- MATLAB source: [`kernel/includes/parallel_profiler_start.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/parallel_profiler_start.m)
- Existing Wiki page: https://spindynamics.org/wiki/index.php?title=parallel_profiler_start.m
