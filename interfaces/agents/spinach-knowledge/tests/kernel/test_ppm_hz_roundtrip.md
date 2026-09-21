# tests/kernel/test_ppm_hz_roundtrip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ppm_hz_roundtrip.m`
- Signature: `result=test_ppm_hz_roundtrip()`
- Total lines: 45

## Purpose

Tests chemical shift and frequency conversion. Syntax: result=test_ppm_hz_roundtrip()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `spin()`, `ppm2hz()`, `test_close()`, `hz2ppm()`.
