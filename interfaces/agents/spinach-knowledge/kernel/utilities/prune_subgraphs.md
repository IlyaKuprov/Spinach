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
