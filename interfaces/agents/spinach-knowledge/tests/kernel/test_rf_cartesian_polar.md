# tests/kernel/test_rf_cartesian_polar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_rf_cartesian_polar.m`
- Signature: `result=test_rf_cartesian_polar()`
- Total lines: 51

## Purpose

Tests RF Cartesian and polar waveform conversion. Syntax: result=test_rf_cartesian_polar()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `polar2cartesian()`, `cartesian2polar()`, `test_close()`.
