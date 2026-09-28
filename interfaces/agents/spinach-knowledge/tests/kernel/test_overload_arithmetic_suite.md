# tests/kernel/test_overload_arithmetic_suite.m

- Signature: `result=test_overload_arithmetic_suite()`

## Purpose

Tests cheap overload arithmetic for cell, struct, RCV, and polyadic classes.

## Physical / mathematical content

## Numerical / algorithmic content

The test compares overload arithmetic on small examples with explicit MATLAB matrix references, using tolerances of `1e-15` for cell, struct, and RCV checks and `1e-14` for most polyadic arithmetic checks.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announce the test target and initialize the regression result.
- Check cell addition, subtraction, scalar addition, elementwise scaling, left and right matrix multiplication, sparse totals, inflation, and complex conversion.
- Check recursive struct addition and left multiplication, including nested fields.
- Compare RCV sparse storage, arithmetic, transposition, multiplication, concatenation, and size against MATLAB sparse references.
- Compare polyadic storage and arithmetic with explicitly opened sums of Kronecker products, including size, nonzero counts, addition, scaling, vector multiplication, transposition, and Kronecker products.