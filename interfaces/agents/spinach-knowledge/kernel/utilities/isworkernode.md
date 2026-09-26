# kernel/utilities/isworkernode.m

- Signature: `answer=isworkernode()`

## Purpose

Returns true if executed inside a parfor or spmd block. This function is used in the internal decision making of Spinach kernel: certain algorithms are switched to their serial ver- sions if the calculation is already running inside some par- allel loop. Syntax: answer=isworkernode()

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
