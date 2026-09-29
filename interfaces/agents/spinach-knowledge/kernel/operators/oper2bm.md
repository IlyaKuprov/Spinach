# kernel/operators/oper2bm.m

- Signature: `[states,coeffs]=oper2bm(A)`
- Direct source: [kernel/operators/oper2bm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/oper2bm.m)
- Wiki: [oper2bm.m](https://spindynamics.org/wiki/index.php?title=oper2bm.m)

## Purpose and basis

Expands a numeric square matrix `A` into the bosonic-monomial operator basis for `n=size(A,1)`. The basis generator [`boson_mono.m`](boson_mono.md) constructs `B(k,q)=(A.c^k)*(A.a^q)` for `k,q=0,...,n-1`, using the creation and annihilation matrices returned by [`weyl.m`](weyl.md). In the struct returned by `weyl`, the ladder entries carry square-root occupation factors; its source records `A.c*A.a=A.n`, `[A.n,A.c]=A.c`, `[A.n,A.a]=-A.a`, and `[A.a,A.c]=A.u` except at the truncated edge, where that commutator's edge element is `1-n`. The monomials obey `[N,B(k,q)]=(k-q)*B(k,q)`. `oper2bm` adds no further scaling to these monomials.

There are `n^2` basis matrices, each `n×n`. The cell array is populated as `BM{k+1,q+1}`; MATLAB column-major linear indexing therefore varies `k` fastest, then `q`. The output index vector is the zero-based linear position `0,...,n^2-1` in that list. Use [`lin2kq.m`](../indexing/lin2kq.md) to convert these Spinach BM indices to `K,Q` indices.

## Coefficients

The routine forms the full Gram matrix `S(i,j)=hdot(BM{i},BM{j})` and overlap vector `v(i)=hdot(BM{i},A)`, then solves `coeffs=S\full(v)`. Here `hdot(X,Y)` is the Frobenius product `sum(sum(conj(X).*Y))` ([`hdot.m`](../utilities/hdot.md)). Thus coefficients are obtained by solving with the computed basis overlaps; the code does not assume unit-norm monomials or divide each overlap by a common norm.

`states` and `coeffs` are returned as corresponding column vectors after retaining only entries satisfying `abs(coeffs)>10*eps('double')`. Each returned state is a basis-list index, not a Hilbert-space state or a propagator index.
