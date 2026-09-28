# kernel/operators/weyl.m

- Signature: `A=weyl(nlevels)`

## Purpose

Construct sparse Weyl boson operators for a bosonic mode truncated to a specified number of population levels.

## Physical / mathematical content

The operators obey `A.c*A.a=A.n`, `[A.n,A.c]=A.c`, `[A.n,A.a]=-A.a`, and `[A.a,A.c]=A.u`, except at the truncation edge state, where the `[A.a,A.c]` element is `1-nlevels`. This exception is unavoidable for finite truncation.

## Numerical / algorithmic content

The operators are constructed as sparse `nlevels`-by-`nlevels` matrices and declared complex at build time to avoid expensive reallocations later.

## Parameters / inputs

- `nlevels` — a positive integer specifying the number of population levels.

## Outputs

- `A.u` — unit operator.
- `A.c` — creation operator.
- `A.a` — annihilation operator.
- `A.n` — population number operator.

## Implementation structure

The function validates `nlevels`, then constructs the creation operator on the lower diagonal, the population number operator on the main diagonal, the annihilation operator on the upper diagonal, and the unit operator as a sparse identity matrix.

<https://spindynamics.org/wiki/index.php?title=weyl.m>