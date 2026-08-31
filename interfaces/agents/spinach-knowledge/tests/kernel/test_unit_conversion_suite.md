# tests/kernel/test_unit_conversion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_unit_conversion_suite.m`
- Signature: `result=test_unit_conversion_suite()`
- Total lines: 62

## Purpose

Tests scalar unit-conversion functions. Syntax: result=test_unit_conversion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Unit-conversion functions\n')`.
- Lines 20-23: State the conversion target of the test; implemented by `result=new_test_result('kernel/unit_conversion_suite', 'Unit-conversion functions', 'unit conversion helpers must implement their defining physical constants.')`.
- Lines 25-26: Check Hartree energy conversion; implemented by `hartree=[0 1 2.5]`.
- Lines 30-31: Check inverse-centimetre and frequency conversion; implemented by `icm=[0 1 12.5]`.
- Lines 38-39: Check electron-field hyperfine conversions; implemented by `hfc_gauss=[0 10 25]`.
- Lines 50-51: Check milliTesla to Hz conversion; implemented by `hfc_mt=[0 1 3.5]`.
- Lines 56-57: Check Lorentzian linewidth conversion; implemented by `fwhm=[1 2.5 10]`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/unit_conversion_suite', 'Unit-conversion functions', 'unit conversion helpers must implement their defining physical constants.')`.
- Lines 26: computes `hartree` using `hartree=[0 1 2.5]`.
- Lines 31: computes `icm` using `icm=[0 1 12.5]`.
- Lines 32: computes `hz` using `hz=100*299792458*icm`.
- Lines 39: computes `hfc_gauss` using `hfc_gauss=[0 10 25]`.
- Lines 40: computes `g` using `g=2.0023193043622`.
- Lines 41: computes `muB` using `muB=9.274009994e-24`.
- Lines 42: computes `hbar` using `hbar=1.054571628e-34`.
- Lines 43: computes `conv` using `conv=1e-10*g*muB/(hbar*2*pi)`.
- Lines 44: computes `hfc_mhz` using `hfc_mhz=conv*hfc_gauss`.
- Lines 51: computes `hfc_mt` using `hfc_mt=[0 1 3.5]`.
- Lines 52: computes `hfc_hz` using `hfc_hz=1e-3*g*muB*hfc_mt/(hbar*2*pi)`.
- Lines 57: computes `fwhm` using `fwhm=[1 2.5 10]`.
- Lines 59: computes `'Lorentzian full width at half maximum corresponds to R2` using `'Lorentzian full width at half maximum corresponds to R2=pi*FWHM')`.

## Outputs

- result -regression test result with explanatory messages
- The test checks physically defined unit conversions used across magnetic
- resonance: Hartree to J/mol, cm^-1 to Hz, field/frequency conversions,
- and Lorentzian linewidth to R2.

## Implementation structure

- Tests scalar unit-conversion functions. Syntax:
- result=test_unit_conversion_suite()
- result -regression test result with explanatory messages
- The test checks physically defined unit conversions used across magnetic
- resonance: Hartree to J/mol, cm^-1 to Hz, field/frequency conversions,
- and Lorentzian linewidth to R2.
- Announce the test target
- State the conversion target of the test
- Check Hartree energy conversion
- Check inverse-centimetre and frequency conversion
- Check electron-field hyperfine conversions
- Check milliTesla to Hz conversion

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `hartree2joule()`, `icm2hz()`, `hz2icm()`, `gauss2mhz()`, `mhz2gauss()`, `mt2hz()`, `fwhm2rlx()`.
