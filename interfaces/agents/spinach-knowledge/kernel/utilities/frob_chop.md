# kernel/utilities/frob_chop.m

- Signature: `r=frob_chop(s,tol)`

## Purpose

Truncates SVD decomposition to the user-specified threshold in the Frobenius norm. Syntax: r=frob_chop(s,tol)

## Physical / mathematical content

- Chooses a retained rank from singular-value vector `s` so the discarded tail has Frobenius norm below the requested tolerance `tol`.

## Numerical / algorithmic content

## Parameters / inputs

- s -a vector of singular values for a matrix,
- in descending order
- tol -truncation threshold

## Outputs

- r -the number of singular values to keep

## Implementation structure

- Converts `s` to a real column vector, zeros numerical noise below `numel(s)*eps*max(abs(s))`, and clamps remaining negative entries to zero. It accumulates squared singular values from the smallest upward, then returns the number retained so the discarded tail remains below `tol`; returns `0` if no singular values need retaining.
