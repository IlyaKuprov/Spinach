# kernel/utilities/prune_subgraphs.m

- Signature: `subgraphs=prune_subgraphs(subgraphs)`

## Purpose

Removes subgraphs that are contained entirely within other subgraphs. Syntax: subgraphs=prune_subgraphs(subgraphs)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
