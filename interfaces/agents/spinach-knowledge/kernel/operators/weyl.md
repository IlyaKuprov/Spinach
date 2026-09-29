# kernel/operators/weyl.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/weyl.m
Wiki: https://spindynamics.org/wiki/index.php?title=weyl.m

## Purpose and basis

`weyl(nlevels)` returns sparse, complex matrices for one truncated bosonic mode. `nlevels` must be a positive real integer. Each output is `nlevels`-by-`nlevels`; the number operator has diagonal entries `0,1,...,nlevels-1`, so matrix position `j` corresponds to population `j-1` in MATLAB one-based indexing.

## Operators and normalisation

The returned fields are `A.u` (identity), `A.c` (creation), `A.a` (annihilation), and `A.n` (number). The source constructs `A.c` from `sqrt(1:nlevels)` on diagonal offset `-1`, `A.a` from `sqrt(0:(nlevels-1))` on offset `+1`, and `A.n` from `0:(nlevels-1)` on the main diagonal; `A.u` is `speye(nlevels)`. This records the literal vectors and offsets supplied to `spdiags`.

The source documents the normalisation relations `A.c*A.a=A.n`, `[A.n,A.c]=A.c`, `[A.n,A.a]=-A.a`, and `[A.a,A.c]=A.u`, with the stated finite-cutoff exception: at the edge state, the `[A.a,A.c]` element is `1-nlevels`. The truncation therefore modifies that commutator at the highest retained population. These finite matrices are not a propagator or evolution generator.
