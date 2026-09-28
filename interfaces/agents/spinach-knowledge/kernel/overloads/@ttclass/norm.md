# kernel/overloads/@ttclass/norm.m

- Signature: `ttnorm=norm(ttrain,norm_type) %#NORMOK`

## Purpose

Computes the norm of the matrix represented by a tensor train. Only the Frobenius norm is available; the 1-norm, infinity norm and 2-norm raise errors.

## Physical / mathematical content

For `norm_type='fro'`, the function packs the train, orthogonalises it with `ttort(ttrain,-1)`, and computes `abs(ttrain.coeff)*norm(ttrain.cores{1,1}(:),2)`.

## Numerical / algorithmic content

The result is the Frobenius norm obtained from the packed, orthogonalised representation and the Euclidean norm of its first core.

## Parameters / inputs

- ttrain - a tensor train representation of a matrix
- norm_type:
  - `1`, `inf`, and `2` are not available for ttclass
  - `'fro'` returns the Frobenius norm

## Outputs

- ttnorm - a nonnegative real number

## Implementation structure

- For `'fro'`, call `pack`, then `ttort` with direction `-1`, and evaluate the scaled 2-norm of the first core.
- Raise an error for the unsupported 1-, infinity-, and 2-norm cases.
