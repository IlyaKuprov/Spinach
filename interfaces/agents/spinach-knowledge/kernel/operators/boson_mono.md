# kernel/operators/boson_mono.m

- Signature: `B=boson_mono(nlevels)`

## Purpose

Construct a cell array of monomials in the creation and annihilation operators for a truncated bosonic mode.

## Physical / mathematical content

The returned operators are `B(k,q)=(Cr^k)*(An^q)` for `k,q=0,...,nlevels-1`, where `Cr` and `An` are the creation and annihilation generators from `weyl(nlevels)`. They satisfy `[N,B(k,q)]=(k-q)*B(k,q)`.

## Numerical / algorithmic content

The function generates all `nlevels^2` pairs and orders them by increasing `k+q`; within each equal-sum group, `k` decreases. For `nlevels=3`, the pair order is `(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(2,1),(1,2),(2,2)`.

## Parameters / inputs

- `nlevels` - positive integer number of bosonic ladder population levels; the indices `k` and `q` run from 0 to `nlevels-1`. The source checks that `nlevels` is numeric, real, scalar, at least 1, and integer-valued.

## Outputs

- `B` - cell array containing the `nlevels^2` bosonic monomial matrices in the ordering described above.

## Implementation structure

The routine obtains the truncated creation and annihilation operators from `weyl(nlevels)`, forms every ordered product `(A.c^k)*(A.a^q)`, then reorders the cells by increasing index sum and decreasing `k` within each sum.

## Reference

- [Spin Dynamics documentation for `boson_mono.m`](https://spindynamics.org/wiki/index.php?title=boson_mono.m)
