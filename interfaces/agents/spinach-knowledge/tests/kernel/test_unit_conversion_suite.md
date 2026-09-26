# tests/kernel/test_unit_conversion_suite.m

- Signature: `result=test_unit_conversion_suite()`

## Purpose

Tests scalar unit-conversion functions. Syntax: result=test_unit_conversion_suite()

## Physical / mathematical content

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
