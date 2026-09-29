# tests/kernel/test_unit_conversion_suite.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_unit_conversion_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_unit_conversion_suite.m)

## Purpose

Regression test for the scalar unit-conversion helper functions used across Spinach. The suite verifies that the conversion helpers implement their defining physical constants, covering:

- Hartree to J/mol conversion (`hartree2joule`),
- inverse centimetre to Hz and back (`icm2hz`, `hz2icm`),
- electron-field hyperfine conversions between Gauss and MHz (`gauss2mhz`, `mhz2gauss`),
- milliTesla to Hz conversion (`mt2hz`),
- Lorentzian full width at half maximum to transverse relaxation rate (`fwhm2rlx`).

## Behaviour

The function announces the test target with `fprintf('TESTING: Unit-conversion functions\n')` and initialises a regression test result object via `new_test_result` with the identifier `kernel/unit_conversion_suite`, the description `Unit-conversion functions`, and the requirement that `unit conversion helpers must implement their defining physical constants.`

Each check is performed with `test_close`, which compares the function output against a reference computed from the defining physical constants, using element-wise relative and absolute tolerances:

- `hartree2joule`: input `[0 1 2.5]` is compared against `2625499.62*hartree` with tolerances `1e-10` and `1e-15`; the message states that one Hartree is 2625499.62 J/mol in the Spinach convention.
- `icm2hz`: input `[0 1 12.5]` is compared against `100*299792458*icm` with tolerances `1e-6` and `1e-15`; the message states that one inverse centimetre is `c*100` Hz.
- `hz2icm`: the frequency vector from the previous check is converted back and compared against the original inverse-centimetre vector with tolerances `1e-12` and `1e-12`; the message states that `hz2icm` is the algebraic inverse of `icm2hz`.
- `gauss2mhz`: input `[0 10 25]` Gauss is compared against `conv*hfc_gauss`, where `conv=1e-10*g*muB/(hbar*2*pi)` with `g=2.0023193043622`, `muB=9.274009994e-24`, and `hbar=1.054571628e-34`; tolerances are `1e-12` and `1e-12`; the message states that Gauss hyperfine units are converted through the electron Zeeman frequency.
- `mhz2gauss`: the MHz values from the previous check are converted back and compared against the original Gauss vector with tolerances `1e-12` and `1e-12`; the message states that `mhz2gauss` is the algebraic inverse of `gauss2mhz`.
- `mt2hz`: input `[0 1 3.5]` milliTesla is compared against `1e-3*g*muB*hfc_mt/(hbar*2*pi)` with tolerances `1e-6` and `1e-12`; the message states that milliTesla hyperfine units are converted to linear frequency through `g*muB/hbar`.
- `fwhm2rlx`: input `[1 2.5 10]` is compared against `pi*fwhm` with tolerances `1e-15` and `1e-15`; the message states that Lorentzian full width at half maximum corresponds to `R2=pi*FWHM`.

The `gauss2mhz`, `mhz2gauss`, and `mt2hz` checks pass the electron g-factor `g` as a second argument to the conversion functions.

## Inputs and outputs

- **Inputs**: none. The function takes no arguments.
- **Outputs**: `result` — regression test result object with explanatory messages, accumulated through successive `test_close` calls.

## References

- [Spinach on GitHub](https://github.com/IlyaKuprov/Spinach)
