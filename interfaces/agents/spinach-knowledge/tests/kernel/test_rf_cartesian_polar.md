# tests/kernel/test_rf_cartesian_polar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_rf_cartesian_polar.m`
- Signature: `result=test_rf_cartesian_polar()`
- Total lines: 51

## Purpose

Tests RF Cartesian and polar waveform conversion. Syntax: result=test_rf_cartesian_polar()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: RF Cartesian/polar conversion\n')`.
- Lines 19-22: State the pulse-control target of the test; implemented by `result=new_test_result('kernel/rf_cartesian_polar', 'RF Cartesian/polar conversion', 'Cartesian and polar RF controls must describe the same complex waveform.')`.
- Lines 24-25: Define a waveform away from the zero-amplitude singularity; implemented by `r=[1.0 2.0 3.0]`.
- Lines 30-31: Convert to Cartesian and back; implemented by `[x,y,Dx,Dy]=polar2cartesian(r,p,Dr,Dp)`.
- Lines 34-36: Check coordinate round-trip; implemented by `result=test_close(result,'RF x coordinate',x,r.*cos(p),1e-15,1e-15, 'x control is amplitude times cosine phase')`.
- Lines 44-46: Check gradient round-trip; implemented by `result=test_close(result,'amplitude gradient round-trip',Dr_back,Dr,1e-12,1e-12, 'gradient components transform by the chain rule')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/rf_cartesian_polar', 'RF Cartesian/polar conversion', 'Cartesian and polar RF controls must describe the same complex waveform.')`.
- Lines 25: computes `r` using `r=[1.0 2.0 3.0]`.
- Lines 26: computes `p` using `p=[0.0 pi/3 -pi/2]`.
- Lines 27: computes `Dr` using `Dr=[4.0 5.0 6.0]`.
- Lines 28: computes `Dp` using `Dp=[0.5 -0.25 0.75]`.
- Lines 31: computes `[x,y,Dx,Dy]` using `[x,y,Dx,Dy]=polar2cartesian(r,p,Dr,Dp)`.
- Lines 32: computes `[r_back,p_back,Dr_back,Dp_back]` using `[r_back,p_back,Dr_back,Dp_back]=cartesian2polar(x,y,Dx,Dy)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that RF amplitude/phase coordinates round-trip to X/Y
- controls and that gradients transform by the chain rule.

## Implementation structure

- Tests RF Cartesian and polar waveform conversion. Syntax:
- result=test_rf_cartesian_polar()
- result -regression test result with explanatory messages
- The test checks that RF amplitude/phase coordinates round-trip to X/Y
- controls and that gradients transform by the chain rule.
- Announce the test target
- State the pulse-control target of the test
- Define a waveform away from the zero-amplitude singularity
- Convert to Cartesian and back
- Check coordinate round-trip
- Check gradient round-trip

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `polar2cartesian()`, `cartesian2polar()`, `test_close()`.
