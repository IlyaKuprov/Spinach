# kernel/overloads/@ttclass/rdivide.m

Signature: `a=rdivide(a,b)`

Source: [kernel/overloads/@ttclass/rdivide.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/rdivide.m) · Wiki: [ttclass/rdivide.m](https://spindynamics.org/wiki/index.php?title=ttclass/rdivide.m)

The method proceeds only when `a` is a `ttclass` and `b` is scalar. Although the source comment describes a numeric scalar, the guard itself checks only `isscalar(b)`; it does not separately require numeric, real, finite, or nonzero `b`. If either check fails, the method raises its stated error; a scalar of another class passes the guard and is subjected to the divisions.

It divides the coefficient vector by `b` and the tolerance vector by `abs(b)`; the core arrays are left unchanged. The tensor train stays in its stored representation rather than being expanded to a full array. For complex `b`, coefficients use `b` itself, not its conjugate; only the tolerance scaling uses its magnitude. There is no zero-divisor check, and no full tensor materialisation occurs in this method.
