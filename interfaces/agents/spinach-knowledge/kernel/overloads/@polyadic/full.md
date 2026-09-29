# kernel/overloads/@polyadic/full.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/full.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/full.m)

`full(p)` materialises the represented matrix. It first recursively converts nested polyadics in cores, prefixes, and suffixes. For each outer core term it starts with the full first factor, appends later factors with `kron` in their stored order, and adds that Kronecker-product term to the accumulated matrix. The row and column counts are products of the factor dimensions in the first term; the source does not check that later terms have matching dimensions.

After summing terms, it multiplies by prefixes in stored product order and then suffixes in stored product order. The accumulator is created with `zeros`, the first factor and final result are explicitly made full, and the source documentation states that full arithmetic is used even when cores are sparse. This is materialisation, not a contraction-only operation; it returns an ordinary full matrix rather than a polyadic. The overload contains no broadcasting or dimension-validation step beyond the matrix and `kron` operations it invokes.

Source comment: [DOI 10.1016/j.evolhumbehav.2017.04.001](http://dx.doi.org/10.1016/j.evolhumbehav.2017.04.001).