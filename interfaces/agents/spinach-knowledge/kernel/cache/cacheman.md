# kernel/cache/cacheman.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/cacheman.m`
- Signature: `cacheman(spin_system) %#NHEAD`
- Total lines: 103

## Purpose

Cache management heuristics. Looks after the scratch folder and prevents it from filling up the disk. Do not call directly. The function inspects the scratch folder and deletes any files that are older than the threshold (default is 365 days) speci- fied in spin_system.tols.cache_mem field.

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Cache management heuristics. Looks after the scratch folder and
- prevents it from filling up the disk. Do not call directly.
- The function inspects the scratch folder and deletes any files
- that are older than the threshold (default is 365 days) speci-
- fied in spin_system.tols.cache_mem field.
- Check consistency
- Get parallel pool directory
- Calculate the time horizon
- Look into the scratch directory
- Delete anything that is out of date
- Report to the user
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gcp()`, `save()`, `delete()`, `dir()`, `any()`, `report()`, `dir_cont()`, `strcmp()`, `rmdir()`, `num2str()`, `isfield()`, `ischar()`, `isscalar()`, `exist()`.
