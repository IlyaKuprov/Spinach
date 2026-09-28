# kernel/utilities/scomponents.m

- Signature: `sci=scomponents(A)`

## Purpose

Finds the strongly connected components of a graph using David Gleich's implementation of Tarjan's algorithm. Reference: http://dx.doi.org/10.1137/0201010

## Syntax

```matlab
sci=scomponents(A)
```

## Parameters / inputs

- `A` — a square logical matrix with 1 for connected nodes in the graph. The function rejects inputs that are not square logical matrices.

## Outputs

- `sci` — a column vector of integers specifying the strongly connected component to which each graph node belongs.

## Implementation structure

The function checks the input with `grumble(A)`, obtains CSR indices with `[rp,ci]=sparse2csr(sparse(A))`, and runs Tarjan's algorithm. Component numbers start at 1 and are assigned as components are found.

Contacts: dgleich@purdue.edu; ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=scomponents.m>