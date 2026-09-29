# kernel/operators/oper2ist.m

- Signature: `[states,coeffs]=oper2ist(A)`
- Direct source: [kernel/operators/oper2ist.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/oper2ist.m)
- Wiki: [oper2ist.m](https://spindynamics.org/wiki/index.php?title=oper2ist.m)

## Purpose and tensor basis

Expands a numeric square matrix `A` into the single-spin irreducible spherical-tensor (IST) basis for multiplicity `n=size(A,1)`. [`irr_sph_ten.m`](irr_sph_ten.md) generates tensors `T(k,m)` satisfying `[Lz,T(k,m)]=m*T(k,m)`, for ranks `k=0,...,n-1` and projections `m=k,...,-k`. Rank zero is the identity. For `k>0`, the top component is initialised as `(-1)^k*2^(-k/2)*L.p^k`, with `L=pauli(n)`; successive lower projections are built from the commutator with `L.m`, divided by `sqrt((k+m)*(k-m+1))` for the previous projection `m`.

There are `n^2` tensors, each `n×n`. In the one-argument generator, list order is increasing rank and, within a rank, decreasing projection. The output `states` are zero-based list positions `0,...,n^2-1`: for example, indices `0,1,2,3` map to `(L,M)=(0,0),(1,1),(1,0),(1,-1)`. Use [`lin2lm.m`](../indexing/lin2lm.md) to convert the returned IST indices to rank and projection.

## Coefficients

For each basis tensor `X`, the implementation computes `hdot(X,A)/hdot(X,X)`, where `hdot(X,Y)=sum(sum(conj(X).*Y))` is the Frobenius product ([`hdot.m`](../utilities/hdot.md)). Each coefficient is therefore normalised by that particular tensor's computed squared Frobenius norm; the code does not assume that every basis tensor has unit norm. It retains entries satisfying `abs(coeffs)>10*eps('double')` and returns their matching indices and coefficients as column vectors.
