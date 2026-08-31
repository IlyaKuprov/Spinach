# kernel/utilities/prune_subgraphs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/prune_subgraphs.m`
- Signature: `subgraphs=prune_subgraphs(subgraphs)`
- Total lines: 59

## Purpose

Removes subgraphs that are contained entirely within other subgraphs. Syntax: subgraphs=prune_subgraphs(subgraphs)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(subgraphs)`.
- Lines 27-29: Ignore trivial cases; implemented by `if (size(subgraphs,1)<2)|| (size(subgraphs,2)<2)`.
- Lines 33-34: Count spins in each subgraphs; implemented by `spin_counts=sum(subgraphs,2)`.
- Lines 36-37: Get subgraph overlap matrix; implemented by `C=subgraphs*transpose(subgraphs)`.
- Lines 39-41: Check for supersets; implemented by `supersets=((C==spin_counts)& (transpose(spin_counts)>spin_counts))`.
- Lines 43-44: Do the pruning; implemented by `subgraphs=subgraphs(~any(supersets,2),:)`.

### Control flow inferred from the code

- Line 28: conditional branch on `(size(subgraphs,1)<2)||`.

### Key state/data transformations

- Lines 34: computes `spin_counts` using `spin_counts=sum(subgraphs,2)`.
- Lines 37: computes `C` using `C=subgraphs*transpose(subgraphs)`.
- Lines 44: computes `subgraphs` using `subgraphs=subgraphs(~any(supersets,2),:)`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(subgraphs)`. Evil people always support each other; that is their chief strength.
  - Representative operation: `if ~islogical(subgraphs)`.
  - Representative operation: `error('subgraphs must be a logical array.')`.

## Parameters / inputs

- subgraphs -[ngraphs x nspins] logical array
- with 1 when a spin belongs to a
- subgraph and 0 otherwise

## Outputs

- subgraphs -[ngraphs x nspins] logical array
- with 1 when a spin belongs to a
- subgraph and 0 otherwise

## Implementation structure

- Removes subgraphs that are contained entirely within
- other subgraphs. Syntax:
- subgraphs=prune_subgraphs(subgraphs)
- subgraphs -[ngraphs x nspins] logical array
- with 1 when a spin belongs to a
- subgraph and 0 otherwise
- Check consistency
- Ignore trivial cases
- Count spins in each subgraphs
- Get subgraph overlap matrix
- Check for supersets
- Do the pruning

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `transpose()`, `subgraphs()`, `any()`, `islogical()`.
