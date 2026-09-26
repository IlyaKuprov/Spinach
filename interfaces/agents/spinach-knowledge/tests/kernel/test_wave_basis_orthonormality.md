# tests/kernel/test_wave_basis_orthonormality.m

- Signature: `result=test_wave_basis_orthonormality()`

## Purpose

Tests waveform basis orthonormality. Syntax: result=test_wave_basis_orthonormality()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks that sine, cosine, and Legendre waveform bases returned
- by Spinach have orthonormal columns, as required by pulse optimisation.

## Implementation structure

- Tests waveform basis orthonormality. Syntax:
- result=test_wave_basis_orthonormality()
- result -regression test result with explanatory messages
- The test checks that sine, cosine, and Legendre waveform bases returned
- by Spinach have orthonormal columns, as required by pulse optimisation.
- Announce the test target
- State the numerical target of the test
- Check all supported basis families
- Build a small waveform basis
- Check orthonormality of columns
