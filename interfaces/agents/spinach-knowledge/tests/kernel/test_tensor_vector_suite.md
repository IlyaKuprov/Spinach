# tests/kernel/test_tensor_vector_suite.m

- Signature: `result=test_tensor_vector_suite()`

## Purpose

Tests tensor, vector, distribution, and relaxation utilities. Syntax: result=test_tensor_vector_suite()

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks Hermite splines, skew-normal density in the normal
- limit, Fokker-Planck vector reshaping helpers, correlation-function
- coefficients, tensor isotope-shift helpers, and small spin-system
- tensor extractors.

## Implementation structure

- Tests tensor, vector, distribution, and relaxation utilities. Syntax:
- result=test_tensor_vector_suite()
- result -regression test result with explanatory messages
- The test checks Hermite splines, skew-normal density in the normal
- limit, Fokker-Planck vector reshaping helpers, correlation-function
- coefficients, tensor isotope-shift helpers, and small spin-system
- tensor extractors.
- Announce the test target
- State the utility target of the test
- Check cubic Hermite interpolation on an exactly representable parabola
- Check skew normal density in the zero-skew normal-distribution limit
- Check phantom-to-Fokker-Planck state embedding
