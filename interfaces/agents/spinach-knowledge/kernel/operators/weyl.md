# kernel/operators/weyl.m

- Signature: `A=weyl(nlevels)`

## Purpose

Weyl boson operators (sparse, see below for normalisa- tion convention) for a bosonic mode with a user-speci- fied population number truncation. Syntax: A=weyl(nlevels)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

## Parameters / inputs

- nlevels -an integer specifying the
- number of population levels

## Outputs

- A.u -unit operator
- A.c -creation operator
- A.a -annihilation operator
- A.n -population number operator
- Note: the matrices are normalised to obey the following
- relations for all energy level counts
- A.c*A.a=A.n, [A.n,A.c]=A.c
- [A.n,A.a]=-A.a, [A.a,A.c]=A.u
- except for the edge state at which [A.a,A.c] ele-
- ment is (1-nlevels), this is unavoidable.
- Note: arrays are declared complex at build time to avoid
- expensive reallocation operations later on.

## Implementation structure

- Weyl boson operators (sparse, see below for normalisa-
- tion convention) for a bosonic mode with a user-speci-
- fied population number truncation. Syntax:
- A=weyl(nlevels)
- nlevels -an integer specifying the
- number of population levels
- A.u -unit operator
- A.c -creation operator
- A.a -annihilation operator
- A.n -population number operator
- Note: the matrices are normalised to obey the following
- relations for all energy level counts
