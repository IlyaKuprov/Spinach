# kernel/operators/irr_sph_ten.m

- Signature: `T=irr_sph_ten(mult)` or `T=irr_sph_ten(mult,k)`
- DIRECT source: [kernel/operators/irr_sph_ten.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/irr_sph_ten.m)
- Wiki: [irr_sph_ten.m](https://spindynamics.org/wiki/index.php?title=irr_sph_ten.m)

## Definition and basis ordering

For a selected rank `k`, the function returns a `2*k+1`-element cell array of `mult`-by-`mult` spin operators, ordered by decreasing projection: `m=k,k-1,...,-k`. The source documents the component relation `Lz*T(k,m)-T(k,m)*Lz = m*T(k,m)`.

With one input, `irr_sph_ten(mult)` recursively gathers all ranks `k=0,...,mult-1`. The rank-`k` block occupies cell positions `k^2+1` through `(k+1)^2`, so the full result has `mult^2` cells, with ranks increasing and projections decreasing within each rank.

## Construction and normalisation

For rank zero the sole tensor is `speye(mult)`, with no division by `sqrt(mult)`. For positive rank, `L=pauli(mult)` supplies the raising and lowering matrices. The highest-projection tensor is initialised exactly as `T{1}=((-1)^k)*(2^(-k/2))*L.p^k`. Subsequent components are generated for `n=2,...,2*k+1` by setting `q=k-n+2` and applying

`T{n}=(L.m*T{n-1}-T{n-1}*L.m)/sqrt((k+q)*(k-q+1))`.

Thus the phase and scale of the top component and each ladder normalisation are explicit in the source; no further normalisation is applied afterward. This routine constructs matrix operators; it does not exponentiate them into propagators.

## Inputs

`mult` must be a finite positive integer. In the two-input form, `k` must be a finite integer satisfying `0<=k<mult`. Calls with any number of inputs other than one or two raise an error. The single-input form returns all of the stated rank blocks.
