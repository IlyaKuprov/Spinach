# kernel/cache/cacheman.m

`cacheman(spin_system)` is an internal scratch-cache cleanup helper and should not be called directly. It uses `spin_system.sys.scratch` as the directory and `spin_system.tols.cache_mem` as the age threshold, which must be a finite, non-negative real scalar. The documented default threshold is 365 days. Because the cutoff is computed as `now - cache_mem`, the tolerance is measured in MATLAB serial-date days.

A write probe creates, saves, and removes a temporary MAT-file in the scratch directory. The function errors if it cannot write there or if the scratch directory is missing. It examines direct entries matching `spinach_*` and treats entries with a modification date earlier than the cutoff as stale. Stale files are deleted; stale directories are removed recursively. Individual deletion failures are caught quietly, and successful file and directory removals are reported.

If a parallel pool already exists, the helper obtains its `Cluster.JobStorageLocation` without starting a pool. It will not recursively remove an expired directory whose full path exactly matches that pool directory.

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/cacheman.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=cacheman.m)
