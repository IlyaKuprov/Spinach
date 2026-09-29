# kernel/overloads/@ttclass/plus.m

Signature: `a=plus(a,b)`

Source: [kernel/overloads/@ttclass/plus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/plus.m) · Wiki: [ttclass/plus.m](https://spindynamics.org/wiki/index.php?title=ttclass/plus.m)

Both inputs must be `ttclass` objects. The method requires equal core counts and equal physical-size arrays; it does not accept a scalar operand. It concatenates the operands' coefficient vectors, core-buffer columns, and tolerance vectors in order. It then removes entries with zero coefficients, applying the same selection to cores and tolerances.

This is a buffered, lazy sum: it retains each nonzero component's existing cores and bond ranks rather than adding core entries or forming/recompressing the represented tensor. If every coefficient is zero, the method replaces the result with `0*unit_like(a)`. No conjugation is performed.
