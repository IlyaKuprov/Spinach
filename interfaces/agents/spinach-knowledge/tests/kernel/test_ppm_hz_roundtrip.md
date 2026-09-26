# tests/kernel/test_ppm_hz_roundtrip.m

- Signature: `result=test_ppm_hz_roundtrip()`

## Purpose

Tests chemical shift and frequency conversion. Syntax: result=test_ppm_hz_roundtrip()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks the physical definition nu=delta*1e-6*gamma*B0/(2*pi)
- and verifies that the inverse conversion preserves the sign of gamma.

## Implementation structure

- Tests chemical shift and frequency conversion. Syntax:
- result=test_ppm_hz_roundtrip()
- result -regression test result with explanatory messages
- The test checks the physical definition nu=delta*1e-6*gamma*B0/(2*pi)
- and verifies that the inverse conversion preserves the sign of gamma.
- Announce the test target
- State the physical target of the test
- Define a field and shifts for positive-gamma and negative-gamma nuclei
- Check explicit frequency formula for proton shifts
- Check that negative magnetogyric ratios retain sign
