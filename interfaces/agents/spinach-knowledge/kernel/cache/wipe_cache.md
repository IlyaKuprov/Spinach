# kernel/cache/wipe_cache.m

- Signature: `wipe_cache(spin_system)`
- Implementation: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/wipe_cache.m>

## Contract

Requests cache cleanup using the scratch-folder settings in a Spinach system object. If called without an argument, the implementation obtains one with `bootstrap('hush')`. It then checks that `spin_system.sys.scratch` is present and names an existing directory, reports the request, sets `spin_system.tols.cache_mem = 0`, and calls `cacheman(spin_system)`. The function has no explicit return value.

## What the cleanup call does

The called cache manager computes its age cutoff as `now - spin_system.tols.cache_mem`; with the zero value set here, the cutoff is the current time. It examines entries matching `spinach_*` in `spin_system.sys.scratch` and removes entries whose modification date is strictly earlier than that cutoff. Old directories are removed recursively; individual removal failures are caught quietly. The manager first tests that it can write a temporary MAT-file in the scratch directory and errors if that test fails. It reports the number of stale files and directories it removed.

Accordingly, `wipe_cache` requests zero-age cleanup through `cacheman`; it does not itself issue a blanket deletion of the scratch directory. The entries considered and the results are those of the cache manager's `spinach_*` scan.

## Input

- `spin_system` — Spinach object whose `sys.scratch` field identifies an existing scratch directory. The no-argument path uses `bootstrap('hush')`.

## Reference

- [Spinach Wiki: `wipe_cache.m`](https://spindynamics.org/wiki/index.php?title=wipe_cache.m)
