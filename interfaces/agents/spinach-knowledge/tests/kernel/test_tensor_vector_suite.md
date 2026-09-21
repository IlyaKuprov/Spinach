# tests/kernel/test_tensor_vector_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_tensor_vector_suite.m`
- Signature: `result=test_tensor_vector_suite()`
- Total lines: 136

## Purpose

Tests tensor, vector, distribution, and relaxation utilities. Syntax: result=test_tensor_vector_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_corr_system()`, `local_tensor_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `herm_spline()`, `test_close()`, `snormpdf()`, `phan2fpl()`, `phantom()`, `fpl2phan()`, `fpl2rho()`, `local_corr_system()`, `corrfun()`, `test_true()`, `isequal()`, `logical()`, `shift_iso()`, `local_tensor_system()`, `get_coupling()`.
