# kernel/utilities/poolsize.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/poolsize.m`
- Signature: `n=poolsize()`
- Total lines: 40

## Purpose

Returns the current parallel pool size. Syntax: n=poolsize()

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Get the pool handle; implemented by `p=gcp('nocreate')`.
- Lines 26-27: Query the pool; implemented by `if isempty(p)`.

### Control flow inferred from the code

- Line 27: conditional branch on `isempty(p)`.

### Key state/data transformations

- Lines 24: computes `p` using `p=gcp('nocreate')`.
- Lines 28: computes `n` using `n=0`.

## Parameters / inputs

- none

## Outputs

- n -number of workers in the current
- parallel pool
- Note: when this function is invoked from inside parfor, spmd,
- or asynchronous parallel job, it returns zero.

## Implementation structure

- Returns the current parallel pool size. Syntax:
- n=poolsize()
- none
- n - number of workers in the current
- parallel pool
- Note: when this function is invoked from inside parfor, spmd,
- or asynchronous parallel job, it returns zero.
- Get the pool handle
- Query the pool
- Run a man through a line of fire, and he turns into a seasoned
- wolf; the weak, and in really tough cases unnecessary, intellect
- is replaced by the wise animal instinct.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gcp()`.
