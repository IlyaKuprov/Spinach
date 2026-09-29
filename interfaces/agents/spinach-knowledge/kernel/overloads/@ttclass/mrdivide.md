# kernel/overloads/@ttclass/mrdivide.m

Direct mapped source: [kernel/overloads/@ttclass/mrdivide.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/mrdivide.m) · [existing Spinach Wiki entry](https://spindynamics.org/wiki/index.php?title=ttclass/mrdivide.m)

## Behaviour

For `a/b`, the guard requires `a` to be `ttclass` and `b` to satisfy `isscalar`; it does not additionally constrain `b` to `double`. The operation divides every train coefficient by `b` and every stored tolerance by `abs(b)`. It leaves cores, mode sizes, and ranks unchanged, and neither compresses nor materialises a dense result.

The scalar is not conjugated; `abs(b)` is used only in the tolerance update. There is no explicit zero-divisor guard, so zero follows MATLAB's division behaviour. Other operand combinations raise the method's error.

No numeric example or DOI is present in the source or either existing page.
