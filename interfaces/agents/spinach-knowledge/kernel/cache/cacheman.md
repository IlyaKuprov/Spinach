# kernel/cache/cacheman.m

- Signature: `cacheman(spin_system) %#NHEAD`

## Purpose

Cache management heuristics. Looks after the scratch folder and prevents it from filling up the disk. Do not call directly. The function inspects the scratch folder and deletes any files that are older than the threshold (default is 365 days) speci- fied in spin_system.tols.cache_mem field.

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
