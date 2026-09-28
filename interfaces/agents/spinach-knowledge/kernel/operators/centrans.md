# kernel/operators/centrans.m

- Signature: `A=centrans(mult,type)`

## Purpose

Construct a sparse central-transition spin operator of dimension `mult` and return it as a complex matrix.

## Physical / mathematical content

The operator has nonzero entries only between the two central spin levels, at 1-based indices `mult/2` and `mult/2+1`. The `type` selects the central-transition component: `x`, `y`, `z`, `+` (raising), or `-` (lowering).

For `x`, the two off-diagonal entries are `0.5`; for `y`, the upper and lower entries are `-0.5i` and `+0.5i`; for `z`, the two diagonal entries are `+0.5` and `-0.5`. Type `+` has a single upper off-diagonal entry `1`, and type `-` a single lower off-diagonal entry `1`. All unspecified entries are zero.

## Numerical / algorithmic content

The routine initializes a sparse `mult`-by-`mult` matrix, fills the entries selected by `type`, and converts the result to complex.

## Parameters / inputs

- `mult` - even integer spin multiplicity, at least 2; it sets the matrix dimension. The source checks that it is numeric, real, scalar, at least 2, and even.
- `type` - character operator selector: `x`, `y`, `z`, `+` (raising), or `-` (lowering).

## Outputs

- `A` - sparse complex central-transition operator matrix of size `mult`-by-`mult`.

## Implementation structure

The function validates the even multiplicity and the character selector, writes the central-transition matrix entries at the two central indices, and returns the sparse matrix as complex.

## Reference

- [Spin Dynamics documentation for `centrans.m`](https://spindynamics.org/wiki/index.php?title=centrans.m)
