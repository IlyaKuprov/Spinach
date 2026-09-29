# kernel/cache/ist_product_table.m

`[product_table_left,product_table_right]=ist_product_table(mult)` returns left- and right-action structure coefficients for the irreducible spherical tensor basis `T` built by `irr_sph_ten(mult)`. The input is an integer multiplicity `mult >= 2`; the basis has `mult^2` elements, so each output has dimensions `mult^2`-by-`mult^2`-by-`mult^2`.

The coefficient ordering is:

- `product_table_left(n,m,k)` is the coefficient of `T{k}` in `T{n}*T{m}`.
- `product_table_right(n,m,k)` is the coefficient of `T{k}` in `T{m}*T{n}`.

The implementation computes norms with `hdot` and projects the normalised operator products back onto `T{k}`. This is the left/right multiplicative-action convention associated with Eq. 7.18 in the first edition of IK's book; the source notes that the book omits its normalisation. `lm2lin` and `lin2lm` translate between single and double indices.

For `mult == 2`, the four-by-four-by-four tables are hard-coded and returned without a disk-cache lookup. Other multiplicities use `ist_product_table_<mult>.mat` beside the function: an existing cache is loaded, or tables are built and saved as v7.3 when possible. A failed save raises a warning but leaves the computed outputs available.

[MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/ist_product_table.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ist_product_table.m)
