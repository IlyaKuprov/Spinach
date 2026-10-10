# kernel/overloads/@polyadic/nnz.m

Source: [kernel/overloads/@polyadic/nnz.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/nnz.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/nnz.m)

A `polyadic` stores additive terms as unopened Kronecker products: each `cores{n}` is one term and each `cores{n}{k}` is a matrix factor. `prefix` and `suffix` are separately stored matrix actions around those terms.

## Behaviour

`nnz(p)` initialises a scalar count to zero, then adds `nnz` of every core factor, every prefix matrix, and every suffix matrix. It traverses the stored factors directly; it does not call `full`, build Kronecker products, or count the nonzero entries of the expanded operator. In particular, this is the sum of stored-factor nonzero counts, not a calculation of cancellations or overlap among expanded terms.

The method has no explicit input validation or dimension checks; it expects the `p.cores`, `p.prefix`, and `p.suffix` fields. The scalar count is also used by `plus` to detect a zero-valued stored operand before combining terms.

## Input and output

- `p`: a polyadic object.
- `answer`: scalar sum of the nonzero counts in the stored factors and buffered matrices.
