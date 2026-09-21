# kernel/utilities/scomponents.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/scomponents.m`
- Signature: `sci=scomponents(A)`
- Total lines: 88

## Purpose

Strongly connected components of a graph, David Gleich's imple- mentation of Tarjan's algorithm:

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sparse2csr()`, `root()`, `sci()`, `islogical()`, `ismatrix()`.
