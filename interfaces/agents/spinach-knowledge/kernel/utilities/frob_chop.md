# kernel/utilities/frob_chop.m

## Purpose

Truncates an SVD decomposition to a user-specified tolerance in the Frobenius norm by returning the number of singular values to keep.

## Behaviour

- Syntax: `r=frob_chop(s,tol)`.
- The function first validates its inputs via an internal consistency check (`grumble`).
- Singular values are reshaped into a real column vector; values with magnitude below `numel(s)*eps*max(abs(s))` — the standard numerical rank threshold — are set to zero because the SVD that produced them does not resolve them. Above that threshold the requested tolerance is honoured exactly.
- Negative values are clipped to zero via `s=max(s,0)`.
- The cutting point is found by computing the cumulative sum of squared singular values from the smallest upward (`cumsum(s(end:-1:1).^2)`) and locating the first index where it reaches `tol^2`.
- If no such index exists, the returned rank is 0; otherwise `r=numel(s)-k+1`.

## Inputs and outputs

Inputs:

- `s` — a vector of singular values for a matrix, in descending order; must be a vector of non-negative real numbers (small imaginary parts up to `1e-10*max(abs(s))` and small negative real parts down to `-1e-10*max(abs(s))` are tolerated as "morally equal" to valid values).
- `tol` — truncation threshold; must be a non-negative real scalar.

Outputs:

- `r` — the number of singular values to keep.

Errors are raised with the messages `'tol must be a non-negative real scalar.'` and `'s must be a vector of non-negative real numbers.'`.

## References

- Source: [frob_chop.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/frob_chop.m)
- [Spinach Wiki: frob_chop.m](https://spindynamics.org/wiki/index.php?title=frob_chop.m)
