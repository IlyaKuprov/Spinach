# kernel/utilities/scomponents.m

- Signature: `sci=scomponents(A)`

## Purpose

Strongly connected components of a graph, David Gleich's imple- mentation of Tarjan's algorithm:

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Syntax

```matlab
sci=scomponents(A)
```

## Parameters / inputs

- A -a logical square matrix with 1 for the
- connected nodes in the graph

## Outputs

- sci -a column vector with integers that spe-
- cify the strongly conected component
- that each node of the graph belongs to

## Implementation structure

- Strongly connected components of a graph, David Gleich's imple-
- mentation of Tarjan's algorithm:
- sci=scomponents(A)
- A - a logical square matrix with 1 for the
- connected nodes in the graph
- sci - a column vector with integers that spe-
- cify the strongly conected component
- that each node of the graph belongs to
- Check consistency
- Get the CSR indices
- Run Tarjan's algorithm
- Consistency enforcement
