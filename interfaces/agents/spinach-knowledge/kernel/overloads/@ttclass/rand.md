# kernel/overloads/@ttclass/rand.m

Signature: `tt=rand(tt,ttrank)`

Source: [kernel/overloads/@ttclass/rand.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/rand.m) · Wiki: [ttclass/rand.m](https://spindynamics.org/wiki/index.php?title=ttclass/rand.m)

The input must be a `ttclass` object; `ttrank` must be a numeric, real, scalar positive integer. The method takes the two physical dimensions of each mode from `sizes(tt)`, replaces the core buffer with one train, and fills its cores immediately with MATLAB `rand` values.

For more than one core, core `k` has shape `[r_left, size_k(1), size_k(2), r_right]`: the first left and last right ranks are 1, and every internal bond rank is `ttrank`. A one-core train has shape `[1, size_1(1), size_1(2), 1]`, so `ttrank` does not change its rank. The result has coefficient 1 and tolerance 0. The method generates stored core values, not a lazy random expression; it applies no conjugation.
