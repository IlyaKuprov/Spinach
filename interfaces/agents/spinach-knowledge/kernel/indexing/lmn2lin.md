# kernel/indexing/lmn2lin.m

- Signature: `I=lmn2lin(L,M,N)`
- Direct MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/lmn2lin.m>
- Spin Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=lmn2lin.m>

## Purpose

Maps the rank `L` and left/right indices `M,N` of Wigner D functions to one-based linear indices. Ranks increase; within a rank, `M` decreases from `L` to `-L`, and for each `M`, `N` decreases from `L` to `-L`. Thus (0,0,0) maps to 1, (1,1,1) to 2, and (1,1,0) to 3.

The routine indexes Wigner D labels; it does not calculate Wigner D values or alter a state, operator, Hamiltonian, or interaction. It is an ordering conversion, not a choice of dynamics formalism or unit convention.

## Mapping and guards

Each rank contributes `(2*L+1)^2` entries. The code evaluates `I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1` element-wise. Within rank `L`, the offset follows the descending `M`, then descending `N` order. Inputs must be real numeric integer-valued arrays; the guards require `L>=0`, `abs(M)<=L`, `abs(N)<=L`, and identical sizes for all three arrays. The range checks precede the explicit size-consistency check in the source. The output has the arithmetic array shape.

## Syntax and arguments

`I=lmn2lin(L,M,N)`

- `L` - non-negative integer Wigner-function rank array.
- `M`, `N` - integer indices satisfying `abs(M)<=L` and `abs(N)<=L`.
- `I` - one-based linear indices; `I=1` corresponds to `L=0, M=0, N=0`.
