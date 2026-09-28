# kernel/operators/boson_ortho.m

- Signature: `B=boson_ortho(nlevels)`

## Purpose

Construct an orthogonal set of truncated bosonic operator monomials.

## Physical / mathematical content

The function orthogonalizes the monomials returned by `boson_mono(nlevels)` with respect to the `hdot` inner product. It uses Gram–Schmidt subtraction and does not normalize the resulting operators.

## Numerical / algorithmic content

Each monomial is successively made orthogonal to the preceding operators in the set. The result retains the input set’s number of operators; orthogonality is with respect to `hdot`, not a claim that the operators have unit norm.

## Parameters / inputs

- `nlevels` - positive integer number of bosonic ladder population levels. The source requires a numeric, real, scalar, positive integer.

## Outputs

- `B` - cell array of orthogonalized bosonic operator matrices corresponding to the monomials from `boson_mono(nlevels)`.

## Implementation structure

The routine validates `nlevels`, calls `boson_mono(nlevels)`, and applies an unnormalized Gram–Schmidt procedure using `hdot` to remove components along previously processed operators.

## Reference

- [Spin Dynamics documentation for `boson_ortho.m`](https://spindynamics.org/wiki/index.php?title=boson_ortho.m)
