# tests/kernel/test_unit_conversion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_unit_conversion_suite.m`
- Signature: `result=test_unit_conversion_suite()`
- Total lines: 62

## Purpose

Tests scalar unit-conversion functions. Syntax: result=test_unit_conversion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
