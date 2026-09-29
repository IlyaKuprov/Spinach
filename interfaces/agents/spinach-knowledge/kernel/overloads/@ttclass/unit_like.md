# kernel/overloads/@ttclass/unit_like.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/unit_like.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/unit_like.m)

## Signature

`A=unit_like(A)`

## Purpose and behaviour

Returns an identity in the same representation family as a square matrix or square-operator tensor train.

For a `ttclass` input, the function checks that every core has equal row and column mode sizes. It then creates one local `eye(mode_size)` matrix per core and calls `ttclass(1,core,0)`. The result is a single rank-one tensor train with the same local mode dimensions, coefficient 1, and tolerance 0; input coefficients, buffered columns, and bond ranks are not copied.

For a square sparse matrix it returns `speye(size(A))`; for a square full matrix it returns `eye(size(A))`. A nonsquare matrix or a tensor train with any unequal local row and column mode size raises an error.

## Input and output

- `A` — a full or sparse square matrix, or a `ttclass` representation whose local row and column mode sizes match.
- `A` — an identity matrix or tensor-train identity in the corresponding representation family.
