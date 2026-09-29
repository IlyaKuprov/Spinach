# kernel/overloads/@opium/kron.m

- Signature: `c=kron(a,b)`

## Meaning and behaviour

An `opium` object represents `coeff*I_dim`. When both operands are `opium` objects, the overload returns `opium(a.dim*b.dim,a.coeff*b.coeff)`: it composes the identity factors in the compact representation and does not expand either object into a matrix. This branch performs no explicit dimension-compatibility check.

When just `a` is `opium`, the method calls MATLAB `kron` on `a.coeff*speye(a.dim)` and `b`; when just `b` is `opium`, it calls MATLAB `kron` on `a` and `b.coeff*speye(b.dim)`. Thus the one-object branches explicitly form that scaled identity as a sparse matrix and do not convert the resulting product back to an `opium` object. The product's resulting storage is left to MATLAB's `kron` and the other operand; this source does not promise a dense or sparse result for every operand combination.

The explicit branches test for `opium` operands; the fallback raises an error for operands that reach it. Numeric/numeric products are not implemented by this overload. No cell-operand or broadcasting branch appears here.

## Inputs and output

- Inputs: `a` and `b`, intended as numeric matrices and/or `opium` objects.
- Output: compact `opium` when both operands are `opium`; otherwise the result returned by MATLAB `kron` for the sparse-identity expansion branch.

## Source links

- MATLAB source: [kernel/overloads/@opium/kron.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/kron.m)
- Existing Wiki page: [opium/kron.m](https://spindynamics.org/wiki/index.php?title=opium/kron.m)
