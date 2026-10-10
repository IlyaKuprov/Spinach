# kernel/operators/pauli.m

- Source: [kernel/operators/pauli.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/pauli.m)
- Wiki: [pauli.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=pauli.m)
- Signature: `S=pauli(mult)`

## Purpose

Constructs the sparse spin-operator matrices for one finite spin with Hilbert-space multiplicity `mult`. This is a spin representation, not a bosonic-mode operator constructor.

## Basis and operator definitions

Let `s=(mult-1)/2`. The matrix basis is ordered by magnetic projection `m=s,s-1,...,-s`; all matrices are `mult`-by-`mult`. `S.u` is the identity and `S.z` is diagonal with those projections. `S.p` is the raising operator on the first off-diagonal and `S.m` is its lowering counterpart on the opposite off-diagonal; their ladder entries use the square-root factors `sqrt(s*(s+1)-m*(m+1))` and `sqrt(s*(s+1)-m*(m-1))`, respectively.

The transverse operators are defined as `S.x=(S.p+S.m)/2` and `S.y=(S.p-S.m)/(2i)`; equivalently, the source comments define `S.p=S.x+1i*S.y` and `S.m=S.x-1i*S.y`. The resulting spin matrices satisfy the cyclic commutation relations `[S.x,S.y]=1i*S.z`, `[S.y,S.z]=1i*S.x`, and `[S.z,S.x]=1i*S.y`.

## Construction and inputs

The multiplicity must be a positive real integer. Multiplicities 2 and 3 use explicit spin-half and spin-one matrices; other multiplicities use the general ladder construction above. The returned matrices are sparse and are made complex at construction.

## Output

`S` is a structure containing `u`, `p`, `m`, `x`, `y`, and `z`, each a `mult`-by-`mult` spin operator. This function returns generators/operators, not a time propagator.
