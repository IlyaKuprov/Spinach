# kernel/utilities/swizzle.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/swizzle.m`
- Signature: `tuples=swizzle(index_arrays)`
- Total lines: 58

## Purpose

Flattens out nested index lists and outputs them as an array of tuples in random order. This is useful for flattening nes- ted loops for parallel processing. Syntax: tuples=swizzle(index_arrays)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(index_arrays)`.
- Lines 25-26: Kronecker up the arrays; implemented by `tuples=index_arrays{1}(:)`.
- Lines 32-33: Randomise the tuple list; implemented by `ntuples=size(tuples,1); tuples=tuples(randperm(ntuples),:)`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=2:numel(index_arrays)`.

### Key state/data transformations

- Lines 26: computes `tuples` using `tuples=index_arrays{1}(:)`.
- Lines 33: computes `ntuples` using `ntuples=size(tuples,1); tuples=tuples(randperm(ntuples),:)`.

### Local helper functions

- Line 38: `grumble()` — `function grumble(index_arrays)`.
  - Representative operation: `if ~iscell(index_arrays)`.
  - Representative operation: `error('index_arrays must be a cell array of row vectors.')`.

## Parameters / inputs

- index_arrays -a cell array of row vectors

## Outputs

- tuples -a matrix of tuples in random or-
- der, with tuples listed as rows

## Implementation structure

- Flattens out nested index lists and outputs them as an array
- of tuples in random order. This is useful for flattening nes-
- ted loops for parallel processing. Syntax:
- tuples=swizzle(index_arrays)
- index_arrays -a cell array of row vectors
- tuples -a matrix of tuples in random or-
- der, with tuples listed as rows
- Check consistency
- Kronecker up the arrays
- Randomise the tuple list
- Consistency enforcement
- Acording to a conference rumour, before appointing IK to a tenured

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `tuples()`, `randperm()`, `iscell()`, `isrow()`, `any()`.
