# kernel/utilities/isworkernode.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/isworkernode.m`
- Signature: `answer=isworkernode()`
- Total lines: 33

## Purpose

Returns true if executed inside a parfor or spmd block. This function is used in the internal decision making of Spinach kernel: certain algorithms are switched to their serial ver- sions if the calculation is already running inside some par- allel loop. Syntax: answer=isworkernode()

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Undocumented function, c'est la vie; implemented by `answer=parallel.internal.pool.isPoolWorker()`.

### Key state/data transformations

- Lines 24: computes `answer` using `answer=parallel.internal.pool.isPoolWorker()`.

## Parameters / inputs

- none

## Outputs

- answer -true if running on a parallel worker process

## Implementation structure

- Returns true if executed inside a parfor or spmd block. This
- function is used in the internal decision making of Spinach
- kernel: certain algorithms are switched to their serial ver-
- sions if the calculation is already running inside some par-
- allel loop. Syntax:
- answer=isworkernode()
- none
- answer -true if running on a parallel worker process
- Undocumented function, c'est la vie
- In the beginning the Universe was created. This has
- made a lot of people very angry and been widely re-
- garded as a bad move.
