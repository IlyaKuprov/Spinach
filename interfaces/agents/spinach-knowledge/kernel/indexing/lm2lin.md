# kernel/indexing/lm2lin.m

- Signature: `I=lm2lin(L,M)`
- Direct MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/lm2lin.m>
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=lm2lin.m>

## Purpose

Maps total-angular-momentum rank `L` and projection `M` to a zero-based linear index for spin-state labels. Indices are ordered by increasing rank, and within each rank by decreasing projection. Thus (0,0) maps to 0, (1,1) to 1, and (1,0) to 2.

This is an indexing conversion only: it does not construct or change a spin state, operator, Hamiltonian, or interaction. Its scope is the (L,M) labelling convention; it does not impose a particular Hamiltonian formalism or physical units.

## Mapping and guards

For each element, `I=L.^2+L-M`. The rank-`L` block starts at `L^2`; values for `M=L,L-1,...,-L` occupy consecutive indices through `(L+1)^2-1`. The implementation validates real numeric integer-valued inputs, then requires `L>=0`, `abs(M)<=L`, and equal array sizes. It applies the formula element-wise and preserves the input array shape.

## Syntax and arguments

`I=lm2lin(L,M)`

- `L` - non-negative integer rank array.
- `M` - integer projection array satisfying `abs(M)<=L`, with the same size as `L`.
- `I` - zero-based linear indices; `I=0` corresponds to `L=0, M=0`.
