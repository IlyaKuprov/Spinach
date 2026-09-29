# kernel/includes/parallel_profiler_report.m

- Signature: script include; intended immediately after a profiled `parfor` or `spmd` stage started with `parallel_profiler_start`

## Summary report and worker-data mapping

The summary runs only when `~isworkernode`. It calls `tocBytes(gcp)`, averages the returned rows over workers with `mean(...,1)`, divides both columns by `2^20`, and reports the first value as the average MB received by a worker and the second as the average MB sent back. It also calls `toc()` and reports the elapsed stage time in seconds. The MB labels and direction are the source's own interpretation of those two columns.

## Detailed profiler file

A second guard requires both `~isworkernode` and `'dafuq'` in `spin_system.sys.enable`. Under that guard the include drains `parProfiler.drainLog()` into `parpool_history`, takes the caller name from `dbstack` entry `a(end-1).name`, and builds `filename=[spin_system.sys.scratch filesep datestr(clock,30) '_' a(end-1).name '.mat']`. It saves `parpool_history` to that MAT file. It then calls `drawnow` and reports the saved filename. When either guard is false, this detailed drain-and-save branch does not run.

## Source comment attribution

The source includes a poem attributed to Philip Larkin.

## Source links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/includes/parallel_profiler_report.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=parallel_profiler_report.m)
