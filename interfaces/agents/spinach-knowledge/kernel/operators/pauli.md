# kernel/operators/pauli.m

- Signature: `S=pauli(mult)`

## Purpose

Constructs sparse spin operators for a spin whose Hilbert-space multiplicity is `mult`.

## Physical / mathematical content

The returned operators satisfy `[S.x,S.y]=1i*S.z`, `[S.y,S.z]=1i*S.x`, and `[S.z,S.x]=1i*S.y` for every supported multiplicity. The raising and lowering operators are `S.p=S.x+1i*S.y` and `S.m=S.x-1i*S.y`.

## Numerical / algorithmic content

The matrices are sparse and declared complex at construction. Multiplicities 2 and 3 use hard-coded spin-half and spin-one matrices; other multiplicities are generated from the spin quantum number `(mult-1)/2` and its magnetic projections. `S.x` and `S.y` are formed from `S.p` and `S.m`.

## Parameters / inputs

- `mult` — positive real integer specifying the spin multiplicity.

## Outputs

- `S.u` — unit operator.
- `S.p` — raising operator.
- `S.m` — lowering operator.
- `S.x` — `Sx` observable operator.
- `S.y` — `Sy` observable operator.
- `S.z` — `Sz` observable operator.

## Implementation structure

The input is checked for numeric, real, scalar, integer, and positive value. Matrices are selected or generated according to multiplicity, then `S.x` and `S.y` are computed from the raising and lowering operators.

## References

- [pauli.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=pauli.m)
