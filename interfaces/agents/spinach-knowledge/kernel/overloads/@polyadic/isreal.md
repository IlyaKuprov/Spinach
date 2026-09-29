# kernel/overloads/@polyadic/isreal.m

MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/isreal.m>
Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=polyadic/isreal.m>

## Meaning

The overload checks the stored representation directly. It visits every factor in every `p.cores{n}` term, then every element of `p.prefix` and `p.suffix`, calling `isreal` on each. It returns `false` immediately when any checked element is not real; if all checks pass, it returns `true`.

This is a representation-level predicate, not a test that first contracts or materialises the complete matrix. Consequently the function's criterion is that each stored factor passes `isreal`, rather than a separate check of the final product after possible cancellations. Calls on nested values use MATLAB method dispatch; this function itself contains no separate nested-polyadic traversal branch. It introduces no broadcasting or shape rule.

No orientation-specific data or action is handled here.
