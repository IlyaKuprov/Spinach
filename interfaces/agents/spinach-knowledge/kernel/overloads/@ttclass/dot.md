# kernel/overloads/@ttclass/dot.m

- Signature: `c=dot(a,b)`

## Purpose

Combines two tensor-train (TT) matrix representations through `ctranspose(a)*b`. The source comment calls the result an inner product; the implementation returns the TT multiplication result directly and does not convert it to a MATLAB numeric scalar.

## Representation and result

A `ttclass` stores its train cores in a cell array with one column per train term and stores one coefficient per term. A core at mode `k` has rank and physical dimensions `r_k × m_k × n_k × r_(k+1)`; the represented matrix has dimensions `prod(m_k) × prod(n_k)`. Each term is a chain of rank-coupled TT cores, not a polyadic sum of rank-one mode factors; this overload does not describe a polyadic-object API.

The call first applies the conjugate transpose, which swaps each core's row/column mode dimensions and conjugates its entries, then invokes TT matrix multiplication. The multiplication contracts the matching physical mode and combines the left and right bond ranks, then sums the pairwise train terms with their coefficients. For matching input mode sizes `M × N`, the result is a `ttclass` representing an `N × N` matrix `AᴴB`; if the inputs are column-vector TTs, that represented result is `1 × 1`.

## Checks

Both inputs must be `ttclass` objects, have the same number of cores, and have equal mode-size arrays. The check does not require equal bond ranks. The subsequent TT multiplication checks that its contracted mode sizes agree.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/dot.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/dot.m)
