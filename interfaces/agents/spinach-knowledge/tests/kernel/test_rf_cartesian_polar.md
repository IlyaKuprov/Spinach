# tests/kernel/test_rf_cartesian_polar.m

- Signature: `result=test_rf_cartesian_polar()`

## Purpose

Tests RF Cartesian and polar waveform conversion. Syntax: result=test_rf_cartesian_polar()

## Physical / mathematical content

## Numerical / algorithmic content

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
