# tests/kernel/test_wave_basis_orthonormality.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_wave_basis_orthonormality.m`
- Signature: `result=test_wave_basis_orthonormality()`
- Total lines: 37

## Purpose

Tests waveform basis orthonormality. Syntax: result=test_wave_basis_orthonormality()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Waveform basis orthonormality\n')`.
- Lines 19-22: State the numerical target of the test; implemented by `result=new_test_result('kernel/wave_basis_orthonormality', 'Waveform basis orthonormality', 'pulse waveform basis columns must be orthonormal.')`.
- Lines 24-25: Check all supported basis families; implemented by `basis_types={'sine_waves','cosine_waves','legendre'}`.
- Lines 28-29: Build a small waveform basis; implemented by `B=wave_basis(basis_types{n},5,32)`.
- Lines 31-33: Check orthonormality of columns; implemented by `result=test_close(result,[basis_types{n} ' Gram matrix'],B'*B,eye(5),1e-12,1e-12, 'orthonormal columns give independent waveform coefficients')`.

### Control flow inferred from the code

- Line 26: `for` loop over `n=1:numel(basis_types)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/wave_basis_orthonormality', 'Waveform basis orthonormality', 'pulse waveform basis columns must be orthonormal.')`.
- Lines 25: computes `basis_types` using `basis_types={'sine_waves','cosine_waves','legendre'}`.
- Lines 29: computes `B` using `B=wave_basis(basis_types{n},5,32)`.

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
