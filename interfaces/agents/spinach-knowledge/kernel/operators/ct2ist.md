# kernel/operators/ct2ist.m

- Signature: `[states,coeffs]=ct2ist(mult,type)`

## Purpose

Expand a central-transition spin operator into its contributing irreducible spherical tensor (IST) basis states and coefficients.

## Physical / mathematical content

The input multiplicity `mult` specifies the spin dimension; `type` selects the central-transition operator component `x`, `y`, `z`, `+`, or `-`. The source constructs that operator with `centrans(mult,type)` and passes it to `oper2ist` for expansion.

## Numerical / algorithmic content

The routine returns the IST basis indices and their expansion coefficients. Use `lin2lm` to convert the indices to spherical-tensor `L,M` labels. The source delegates the operator construction and expansion to `centrans` and `oper2ist`; it does not provide a separate closed-form expansion formula.

## Parameters / inputs

- `mult` - even integer spin multiplicity, at least 2; the source checks that it is numeric, real, scalar, at least 2, and divisible by 2.
- `type` - character central-transition selector: `x`, `y`, `z`, `+`, or `-`.

## Outputs

- `states` - contributing states in Spinach IST basis indexing; use `lin2lm` to convert to `L,M` indices.
- `coeffs` - coefficients of the corresponding ISTs in the operator expansion.

## Implementation structure

The function validates `mult` and `type`, obtains the central-transition operator through `centrans(mult,type)`, and calls `oper2ist` to return the state indices and coefficients.

## Reference

- [Spin Dynamics documentation for `ct2ist.m`](https://spindynamics.org/wiki/index.php?title=ct2ist.m)
