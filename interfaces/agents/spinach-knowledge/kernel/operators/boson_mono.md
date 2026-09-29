# kernel/operators/boson_mono.m

Direct source: [kernel/operators/boson_mono.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/boson_mono.m)

- Signature: `B=boson_mono(nlevels)`

## Purpose

Construct the ordered monomials of the truncated bosonic creation and annihilation matrices returned by `weyl(nlevels)`. For each `k,q=0,...,nlevels-1`, the cell contains `A.c^k*A.a^q`, with creation powers on the left and annihilation powers on the right. This function adds no further scalar prefactor; the ladder-matrix convention comes from `weyl`. Its documented number-operator relation is `[N,B(k,q)]=(k-q)*B(k,q)`.

## Ordering and dimensions

The returned `B` is an `nlevels^2-by-1` cell array. Each `B{j}` is an `nlevels-by-nlevels` matrix. The cells are ordered by increasing `k+q` and, within each equal-sum group, decreasing `k`. For `nlevels=3` the 1-based cell-number map is: `1:(0,0), 2:(1,0), 3:(0,1), 4:(2,0), 5:(1,1), 6:(0,2), 7:(2,1), 8:(1,2), 9:(2,2)`.

The function accepts a positive integer `nlevels`. The finite truncation and matrix conventions are those of [`weyl.m`](weyl.md).

## Operator action

The output consists of Hilbert-space operator matrices formed by ordinary matrix products. The routine does not construct left-multiplication, right-multiplication, a commutator superoperator, or a propagator.

## Reference

- [Spin Dynamics documentation for `boson_mono.m`](https://spindynamics.org/wiki/index.php?title=boson_mono.m)
