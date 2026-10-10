# kernel/overloads/@ttclass/minus.m

Direct mapped source: [kernel/overloads/@ttclass/minus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/minus.m) · [existing Spinach Wiki entry](https://spindynamics.org/wiki/index.php?title=ttclass/minus.m)

## Behaviour

For `a-b`, the method represents the difference as concatenated tensor-train components; it does not subtract corresponding cores or call `shrink`. Existing cores and ranks stay with their respective components, while the coefficients of `b`'s components are negated. The component list can therefore grow until a later compression.

Both operands must be `ttclass` objects with the same number of cores and identical mode-size arrays (`sizes(a)==sizes(b)`). Their train counts and ranks are not compared. The method appends the coefficients and tolerances, then removes entries whose coefficients are exactly zero together with their cores and tolerances. If all coefficients are zero, it returns `0*unit_like(a)`.

There is no conjugation: only the second operand's coefficients change sign. No numeric example or DOI is present in the source or either existing page.
