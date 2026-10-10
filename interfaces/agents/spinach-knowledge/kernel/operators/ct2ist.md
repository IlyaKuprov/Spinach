# kernel/operators/ct2ist.m

Direct source: [kernel/operators/ct2ist.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/ct2ist.m)

- Signature: `[states,coeffs]=ct2ist(mult,type)`

## Purpose

Expand the central-transition matrix from `centrans(mult,type)` into the Spinach irreducible spherical tensor (IST) basis. The accepted `type` values are `x`, `y`, `z`, `+`, and `-`; `mult` is even and at least 2. The constructed matrix is `mult-by-mult`.

## Indexing, coefficients, and dimensions

`states` and `coeffs` are corresponding column vectors. The IST basis contains `mult^2` matrices; `oper2ist` numbers them from zero in the order returned by `irr_sph_ten(mult)`. In its MATLAB cell array, rank `L` tensors occupy positions `L^2+1` through `(L+1)^2`; `oper2ist` reports these as zero-based state labels `L^2` through `(L+1)^2-1`. The returned subset is in ascending state order. Use [`lin2lm`](../indexing/lin2lm.md) to convert a state index to spherical-tensor `L,M` labels.

For each basis matrix `T_s`, `oper2ist` computes the coefficient `hdot(T_s,A)/hdot(T_s,T_s)`, where `A` is the central-transition matrix and `hdot(X,Y)=sum(conj(X).*Y,'all')`. It retains only entries with absolute coefficient strictly greater than `10*eps('double')`; the same logical mask is applied to `states` and `coeffs`. Thus the output is a thresholded basis expansion, not a newly normalised operator.

## Operator action

This routine returns basis indices and scalar expansion coefficients, not an action matrix, left/right superoperator, or propagator.

## Reference

- [Spin Dynamics documentation for `ct2ist.m`](https://spindynamics.org/wiki/index.php?title=ct2ist.m)
- Related entries: [`centrans.m`](centrans.md), [`oper2ist.m`](oper2ist.md), [`irr_sph_ten.m`](irr_sph_ten.md)
