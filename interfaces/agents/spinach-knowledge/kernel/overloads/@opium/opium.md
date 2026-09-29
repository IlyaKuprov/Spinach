# kernel/overloads/@opium/opium.m

- MATLAB implementation: [kernel/overloads/@opium/opium.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/opium.m)

## Purpose

M=opium(dim,coeff) constructs an OPIUM (“Object Pretending It is a Unit Matrix”): a compact object representing the scaled identity coeff*I_dim. It stores dim and coeff; it does not construct a dense matrix at creation. Use it where downstream code understands this representation rather than assuming it is a fully materialised MATLAB matrix.

## Inputs and output

Both constructor arguments are required:

- dim: numeric, real, scalar, positive integer dimension. Otherwise the constructor errors with “dim must be a positive integer scalar.” No unit applies.
- coeff: numeric scalar multiplying the identity. The implementation does not require it to be real or finite; its units are those of the represented matrix, if applicable.
- M: an opium object with public properties dim and coeff.

## Representation methods

- sparse(M) materialises the represented matrix as M.coeff*speye(M.dim).
- nnz(M) returns 0 when coeff==0 and 1 otherwise; the method reports one for any nonzero coefficient; it does not return the count in the expanded identity matrix.
- numel(M) is always 1. isnumeric(M) and ismatrix(M) return true. allfinite(M) tests whether the coefficient is finite, and iseye(M) is true exactly when coeff==1.
- conj and conjugate-transpose conjugate the coefficient. gpuArray and gather transfer the coefficient to or from a GPU array; these methods do not by themselves expand the identity.

The class declares gpuArray as an inferior class. No physical units are assigned by this wrapper; dimensional meaning comes from the caller's operator and coefficient.

Source documentation: [opium.m](https://spindynamics.org/wiki/index.php?title=opium/opium.m).
