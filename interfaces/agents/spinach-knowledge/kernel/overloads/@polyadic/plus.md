# kernel/overloads/@polyadic/plus.m

Source: [kernel/overloads/@polyadic/plus.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/plus.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/plus.m)

A `polyadic` stores each additive term as an unopened Kronecker product: an outer `cores` entry is a term, and its inner cells are its matrix factors. Thus `cores={{A,B},{C}}` represents `A kron B + C`. `prefix` and `suffix` hold left- and right-side matrix actions.

## Behaviour and dimensions

The method first checks operand sizes. If both operands are non-scalar, unequal row or column dimensions raise `operands must represent matrices of the same dimension.`. If either operand is scalar, this mismatch check is skipped; the method does not implement general array broadcasting.

It next tests `nnz`: if `a` has zero stored-factor count it returns `simplify(b)`, and conversely for `b`. Otherwise it represents the sum without opening the Kronecker products:

- Numeric matrix plus polyadic: when the polyadic has no prefix or suffix, append the numeric matrix as a one-factor term; with buffered actions, construct two one-factor terms for the operands.
- Polyadic plus numeric matrix: symmetric handling, testing the polyadic operand's prefix and suffix.
- Two polyadics with no prefix or suffix: concatenate their core-term lists, preserving the sum as separate Kronecker terms. If either has buffered actions, wrap the two operands as separate one-factor terms.

The result is passed to `simplify` immediately. The mapped method does not call `full`, but `simplify` can merge adjacent eligible non-`opium` factors with `kron` and unwrap a lone factor from a buffer-free single-term result; some factor products can therefore be materialised eagerly. The mapped method uses ordinary addition and introduces no complex conjugation or non-scalar broadcasting. Use repeated polyadic addition sparingly: sums are buffered as separate terms, and `simplify` does not expand the complete sum; later operations traverse a growing term list and can be slower.

## Inputs and output

- `a` and `b`: operands accepted through MATLAB dispatch as numeric matrices/scalars and/or polyadic objects; the method's branches access polyadic fields when an operand is a polyadic.
- `c`: simplified sum, with the common matrix dimensions when both operands are non-scalar; simplification may unwrap a trivial one-term, buffer-free polyadic to its underlying matrix.
