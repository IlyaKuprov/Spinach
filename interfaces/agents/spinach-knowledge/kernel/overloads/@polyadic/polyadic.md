# kernel/overloads/@polyadic/polyadic.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/polyadic.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic.m)

## Signature

`p=polyadic(cores)`

## Representation

`cores` is a cell array of cell arrays. Each outer cell holds one summand; its inner cell holds the factors of that summand’s Kronecker product. Thus `cores={{A,B,C},{D,E}}` represents `kron(A,kron(B,C)) + kron(D,E)`. Matrix factors may themselves be polyadic objects. The class also initialises empty `prefix` and `suffix` cell arrays for later matrix-composition factors.

The Kronecker products are stored unopened: the source states that multiplicative actions can be performed without opening them, potentially saving orders of magnitude in CPU time. This is a factorised representation, not an eagerly expanded sum.

## Inputs and checks

- `cores`: must be a cell array, and each outer element must also be a cell array; otherwise the constructor raises the corresponding `cores must be a cell array.` or `elements of cores must also be cell arrays.` error.
- After assigning `cores`, the constructor calls `validate(p)` to validate the constructed object. The constructor itself does not state additional dimension rules.

## Output

- `p`: the polyadic object containing the supplied factorisation.

No complex conjugation, scalar broadcasting, or Kronecker-product expansion is performed by this constructor. Related overloads: [prefix](./prefix.md), [simplify](./simplify.md), and [size](./size.md).