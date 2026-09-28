# kernel/cache/wipe_cache.m

- Signature: `wipe_cache(spin_system)`

## Purpose

Forces a wipe of the Spinach cache folder.

## Physical / mathematical content

## Numerical / algorithmic content

- Sets `spin_system.tols.cache_mem` to zero and calls `cacheman(spin_system)`.

## Parameters / inputs

- `spin_system` — Spinach object with the cache folder location in `spin_system.sys.scratch`. If omitted, `bootstrap('hush')` supplies the default object.

## Output

- Attempts to delete all Spinach-specific files in `spin_system.sys.scratch`; this may fail quietly if file system permissions are insufficient.

## Implementation structure

- Checks that `spin_system.sys.scratch` is specified and that the folder exists, then reports the requested cache wipe before invoking cache management.

<https://spindynamics.org/wiki/index.php?title=wipe_cache.m>
