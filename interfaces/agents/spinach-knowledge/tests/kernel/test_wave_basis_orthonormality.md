# tests/kernel/test_wave_basis_orthonormality.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_wave_basis_orthonormality.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_wave_basis_orthonormality.m)

## Purpose

Regression test that verifies the sine, cosine, and Legendre waveform bases returned by Spinach have orthonormal columns, as required by pulse optimisation.

## Behaviour

- Announces the test target by printing `TESTING: Waveform basis orthonormality`.
- Initialises a regression test result via `new_test_result` with test name `kernel/wave_basis_orthonormality`, description `Waveform basis orthonormality`, and the numerical target statement `pulse waveform basis columns must be orthonormal.`
- Iterates over the basis families `sine_waves`, `cosine_waves`, and `legendre`.
- For each family, builds a waveform basis with `wave_basis(basis_type,5,32)`.
- Checks orthonormality by comparing the Gram matrix `B'*B` against `eye(5)` using `test_close` with absolute and relative tolerances of `1e-12` each, logging the check under a label of the form `<basis_type> Gram matrix` with the message `orthonormal columns give independent waveform coefficients`.

## Inputs and outputs

- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** none (the function takes no arguments).

## References

- `wave_basis` — constructs the waveform bases under test.
- `new_test_result` — initialises the regression test result object.
- `test_close` — performs the numerical closeness comparison against the identity Gram matrix.
