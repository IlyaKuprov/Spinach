# kernel/overloads/@polyadic/kron.m

MATLAB source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/kron.m>
Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=polyadic/kron.m>

## Meaning and accepted operands

This overload packages the ordered factors of `kron(a,b)` into a polyadic value rather than explicitly forming the full numeric Kronecker matrix in this function. Its local guard accepts operands that are either instances of `polyadic` or numeric. Although the error text says “matrices or polyadics,” this file checks `isnumeric`, not `ismatrix`; any further shape constraints are left to the constructor and `simplify`.

## Factor placement and composition

If `a` is polyadic, `b` is numeric, and `a` has no prefixes or suffixes, `b` is appended to each core-factor list of `a`. If `a` is numeric, `b` is polyadic, and `b` has no prefixes or suffixes, `a` is prepended to each core-factor list of `b`. These placements retain the order of the requested Kronecker product. In the other cases—including two polyadic operands, or an operand with affixes—the function creates a nested value from `{{a,b}}`, preserving both operands as ordered factors. It then calls `simplify` on the result.

This method does not call `inflate` or state a broadcasting rule. It represents factorised composition; explicit expansion is handled separately by `inflate`. Any structural changes made by `simplify` are outside the logic defined here. There is no orientation-specific operation in this overload.

## Source acknowledgement

The source carries a DOI link in an authorship acknowledgement, retained here as source provenance rather than as support for the implementation: <http://dx.doi.org/10.1046/j.1365-294x.1998.00309.x>.
