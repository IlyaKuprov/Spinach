# kernel/overloads/@polyadic/mtimes.m

Source: [kernel/overloads/@polyadic/mtimes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/mtimes.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/mtimes.m)

`polyadic` storage represents a sum of terms, each term being an unopened Kronecker product of its `cores` factors. `prefix` matrices act on the left and `suffix` matrices on the right; list order is retained. The overload chooses whether to update factors, buffer an operand, or materialise a numeric product.

## Scalar and matrix operands

- When a numeric scalar is a compatible one-element state (`size(A,2)==1`), right multiplication, whether dense or sparse, follows the ordinary numeric-block product path and returns a numeric column. Otherwise, right scalar multiplication remains buffered operator scaling. Left numeric scalar multiplication remains buffered scaling.
- A numeric scalar on the left scales the smallest matrix factor in each term of `B`; operator scaling on the right reuses this left-scaling route. Each route calls `simplify`. The scalar itself is applied as factor multiplication, not as a conjugating operation. The following `simplify` may merge factors as described below.
- A non-scalar sparse numeric `A` is prepended to `B.prefix`; nonscalar sparse numeric `B` is appended to `A.suffix`. These are stored as deferred matrix actions and then simplified.
- A full numeric `A` on the left is routed recursively as `ctranspose(ctranspose(B)*ctranspose(A))`. This is the conjugate-transpose identity for the ordinary product `A*B`, not a request to conjugate the product. The recursive call has a full numeric right operand, so it uses the full-right branch and materialises a full numeric product.
- With a full numeric right operand `B`, the method applies buffered suffix factors in reverse loop order, converts `B` to full, accumulates each core term with `kronm`, then applies prefix factors and returns a full numeric result. This branch materialises the product rather than preserving it as a polyadic.

The operator-scaling and sparse branches do not directly call `full`. Their `simplify` calls may merge adjacent `opium` identities, remove zero terms, flatten eligible nested polyadics, and unwrap a lone numeric factor from a buffer-free single-term result.

## Polyadic by polyadic

The eager core-by-core path is allowed only when `A.suffix` and `B.prefix` are empty, each operand has one core term, both terms have the same number of factors, and each pair has matching inner dimensions (`size(A_i,2)==size(B_i,1)`). For each paired factor not handled as `opium`, every row and column dimension of both factors must be at most `1024`. Before simplification, the result is one term whose factors are `A_i*B_i`; `A.prefix` and `B.suffix` are carried to it. Its matrix dimensions are those of the ordinary product: rows of `A` by columns of `B`. The subsequent `simplify` may merge adjacent identities or unwrap a trivial result.

If any eager-path condition fails, the method calls `suffix(A,B)` and simplifies, retaining `B` as a right-side buffered action rather than directly expanding all Kronecker terms. Simplification may still merge eligible factors as noted above. Branches distinguish numeric and polyadic operands with `isnumeric`, `isa`, `isscalar`, and `issparse` tests. The function has a final explicit error for unsupported operand types, although malformed combinations may fail earlier in a branch-specific field access or matrix operation. Operator scaling preserves the polyadic operand dimensions; matrix multiplication has the ordinary output dimensions `size(A,1)` by `size(B,2)` when the inner dimensions conform. The result may remain polyadic or be unwrapped/materialised by simplification or the full-right path. No general broadcasting rule is implemented.

A function handle core in either operand disables the eager core-by-core multiplication path, irrespective of size; the composition is buffered instead. A scalar factor is absorbed by the smallest matrix core of each term, never by a handle; if a term has only handles, a scalar numeric factor is appended. Block application continues through `kronm`, which calls the handle on the reshaped block.
