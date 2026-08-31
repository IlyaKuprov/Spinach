# tests/kernel/test_ppm_hz_roundtrip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ppm_hz_roundtrip.m`
- Signature: `result=test_ppm_hz_roundtrip()`
- Total lines: 45

## Purpose

Tests chemical shift and frequency conversion. Syntax: result=test_ppm_hz_roundtrip()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Chemical-shift frequency conversion\n')`.
- Lines 19-22: State the physical target of the test; implemented by `result=new_test_result('kernel/ppm_hz_roundtrip', 'Chemical-shift frequency conversion', 'ppm2hz and hz2ppm must implement the Larmor-frequency definition.')`.
- Lines 24-25: Define a field and shifts for positive-gamma and negative-gamma nuclei; implemented by `B0=14.1`.
- Lines 28-29: Check explicit frequency formula for proton shifts; implemented by `hz_ref=1e-6*ppm*(B0*spin('1H')/(2*pi))`.
- Lines 36-37: Check that negative magnetogyric ratios retain sign; implemented by `hz_ref=1e-6*ppm*(B0*spin('15N')/(2*pi))`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/ppm_hz_roundtrip', 'Chemical-shift frequency conversion', 'ppm2hz and hz2ppm must implement the Larmor-frequency definition.')`.
- Lines 25: computes `B0` using `B0=14.1`.
- Lines 26: computes `ppm` using `ppm=[-2.5 0 3.0 12.0]`.
- Lines 29: computes `hz_ref` using `hz_ref=1e-6*ppm*(B0*spin('1H')/(2*pi))`.
- Lines 30: computes `hz_obs` using `hz_obs=ppm2hz(ppm,B0,'1H')`.

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
