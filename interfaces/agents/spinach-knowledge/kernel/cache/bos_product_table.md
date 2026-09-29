# kernel/cache/bos_product_table.m

Enumerates the multiplication structure coefficients for the orthogonalised bosonic monomials produced by `boson_ortho(nlevels)`. The input `nlevels` is a positive integer number of bosonic ladder population levels.

For `nlevels` levels, the basis contains `nlevels^2` monomials and each output is an `nlevels^2`-by-`nlevels^2`-by-`nlevels^2` array. Indices `n`, `m`, and `k` range over that basis. The coefficients follow the source conventions:

- `product_table_left(n,m,k)` is the coefficient of `B{k}` in `B{n}*B{m}`.
- `product_table_right(n,m,k)` is the coefficient of `B{k}` in `B{m}*B{n}`.

The implementation evaluates these coefficients with `hdot` after scaling each monomial by its norm; for example, the left entry is `norms(n)*hdot(B{k}/norms(k),(B{n}/norms(n))*(B{m}/norms(m)))`. This left/right multiplicative-action convention is the one cited at Eq. 7.18 in the first edition of IK's book; the source notes that the book omits its normalisation.

The arrays are cached beside the function in `bos_product_table_<nlevels>.mat`. An existing file is loaded; otherwise the tables are built and saved as v7.3 when possible. A failed cache save warns that the Spinach directory appears write-protected but does not prevent returning the computed tables.

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/bos_product_table.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=bos_product_table.m)
