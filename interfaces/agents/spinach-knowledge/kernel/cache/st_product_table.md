# kernel/cache/st_product_table.m

- Signature: `[pt_left,pt_right] = st_product_table(nlevels)`
- Implementation: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cache/st_product_table.m>

## Contract

For `nlevels` energy levels, the function returns the coefficients of products of the single-transition operator basis `B = sin_tran(nlevels)`. With every index in `1:nlevels^2`, the coefficients are defined by the source conventions `S{n}*S{m} = ... + pt_left(n,m,k)*S{k} + ...` and `S{m}*S{n} = ... + pt_right(n,m,k)*S{k} + ...`. Thus each output is an `nlevels^2-by-nlevels^2-by-nlevels^2` array; the first two indices name the multiplied operators and the third names the projected basis operator.

The implementation evaluates the two ordered products separately: `hdot(B{k}, B{n}*B{m})` for the left table and `hdot(B{k}, B{m}*B{n})` for the right table. The source identifies this convention with Eq. 7.18 of the first edition of IK's book and notes that its normalisation is missing there; it also points to `kq2lin` and `lin2kq` for translating between single and double indices.

## Input

- `nlevels` must be a numeric, real, positive integer scalar.

## Cache behaviour

The function first looks beside its implementation for `st_product_table_<nlevels>.mat`. If present, it loads `pt_left` and `pt_right`; otherwise it constructs the tables and attempts to save those variables there. A save failure is caught: the computed outputs are still returned, with a warning that the installation may be write-protected. This function's disk-cache filename is keyed by `nlevels`.

## Reference

- [Spinach Wiki: `st_product_table.m`](https://spindynamics.org/wiki/index.php?title=st_product_table.m)
