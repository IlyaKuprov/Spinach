# tests/kernel/test_wave_basis_orthonormality.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_wave_basis_orthonormality.m`
- Signature: `result=test_wave_basis_orthonormality()`
- Total lines: 37

## Purpose

Tests waveform basis orthonormality. Syntax: result=test_wave_basis_orthonormality()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `wave_basis()`, `test_close()`.
