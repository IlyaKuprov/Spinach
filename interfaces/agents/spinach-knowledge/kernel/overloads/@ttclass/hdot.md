# kernel/overloads/@ttclass/hdot.m

- Signature: `c=hdot(a,b)`

## Purpose

Computes a scalar Hermitian Frobenius (Hadamard) dot product of two TT matrix representations. It is not a matrix-valued elementwise Hadamard product.

## Core contraction and storage

A `ttclass` can store multiple train terms in the columns of its core-cell array, with one coefficient per term. Each term is a chain of rank-coupled TT cores, rather than a polyadic sum of rank-one mode factors; this overload does not describe a polyadic-object API. A mode-`k` core has dimensions `r_k × m_k × n_k × r_(k+1)`. For each pair of terms, the function initialises the contraction with `conj(a.coeff(na))*b.coeff(nb)`. At each mode it reshapes the `b` core across its left rank, row/column physical dimensions, and right rank; it pairs this with the correspondingly reshaped `a` core and updates the contraction by `core_a'*core_b`. Thus the intermediate is a left-rank-by-right-rank matrix, and the terminal unit ranks reduce it to a scalar. All term-pair scalars are accumulated into `c`.

## Checks and result

Both arguments must be `ttclass` objects, have the same number of cores, and have equal row/column mode sizes at every core. The implementation allows different bond ranks and contracts them pairwise; it does not require identical internal ranks. The result is a numeric scalar, which may be complex for complex-valued trains.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/hdot.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/hdot.m)
