# tests/kernel/test_lowlevel_utilities_suite.m

- Signature: `result=test_lowlevel_utilities_suite()`

## Purpose

Tests cheap deterministic low-level utility functions against explicit answers. Syntax: `result=test_lowlevel_utilities_suite()`.

## Physical / mathematical content

- Checks a Gaussian with unit standard deviation against the standard normal density, and the real and imaginary parts of a zero-phase Lorentzian against `1/(1+x^2)` and `x/(1+x^2)` at `x=[0 1]` for the chosen parameters.
- Checks rank-two rotational spectral density: at zero frequency, `spden(2,10,0)=1/300`; when `omega*tau_c=1`, `spden(2,10,60)=1/600`.

## Numerical / algorithmic content

- Verifies the commutator `AB-BA` and left-to-right nesting by `rocomm`.
- Verifies trace removal and extraction of the part commuting with an operator in its eigenbasis, including degenerate blocks and eigenvalue splittings across large offsets or spectral spans.
- Checks the Frobenius inner product `trace(D'*E)`, anti-diagonal transpose, and helpers that zero selected rows, columns, or a diagonal band.
- Checks reconstruction from leading singular values with `keep_rank` and rank selection by discarded-tail tolerance with `frob_chop`.
- Checks signed and unsigned minimum-integer-type selection at the `127/128` and `255/256` promotion boundaries.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

Creates a test result for `kernel/lowlevel_utilities_suite`, then uses `test_close` for numerical and matrix comparisons and `test_true` for integer-type checks.