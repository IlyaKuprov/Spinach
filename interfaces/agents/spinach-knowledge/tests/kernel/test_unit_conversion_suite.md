# tests/kernel/test_unit_conversion_suite.m

- Signature: `result=test_unit_conversion_suite()`

## Purpose

Tests scalar unit-conversion functions used in magnetic resonance.

## Physical / mathematical content

The suite checks Hartree to J/mol, inverse-centimetre to Hz, field/frequency, and Lorentzian linewidth to R2 conversions. It uses `1 Hartree = 2625499.62 J/mol`, `1 cm^-1 = 100*c Hz`, and `R2 = pi*FWHM`.

## Numerical / algorithmic content

The tests compare each conversion with values calculated from its defining constants using `test_close`, and check that `hz2icm` and `mhz2gauss` invert their corresponding forward conversions.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announce the unit-conversion test and initialize its result.
- Check Hartree energy conversion with inputs `[0 1 2.5]`.
- Check inverse-centimetre to Hz conversion and its inverse with inputs `[0 1 12.5]`.
- Check Gauss to MHz conversion and its inverse with inputs `[0 10 25]`, using the electron g-factor, Bohr magneton, and reduced Planck constant.
- Check milliTesla to Hz conversion with inputs `[0 1 3.5]`.
- Check Lorentzian full width at half maximum to R2 conversion with inputs `[1 2.5 10]`.