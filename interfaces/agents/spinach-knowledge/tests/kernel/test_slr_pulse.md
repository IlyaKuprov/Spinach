# tests/kernel/test_slr_pulse.m

- Signature: `result=test_slr_pulse()`

## Purpose

Tests Shinnar-Le Roux selective excitation pulse design. Syntax: result=test_slr_pulse()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- result -regression test result with explanatory messages
- The test checks waveform units and shape, independent two-level
- propagation, excitation profile selectivity, production-path shaped
- pulse propagation, and representative input validation failures.

## Implementation structure

- Tests Shinnar-Le Roux selective excitation pulse design. Syntax:
- result=test_slr_pulse()
- result -regression test result with explanatory messages
- The test checks waveform units and shape, independent two-level
- propagation, excitation profile selectivity, production-path shaped
- pulse propagation, and representative input validation failures.
- Announce the test target
- State the selective pulse target of the test
- Define a representative selective excitation design
- Generate the production waveform
- Check the output dimensions and finiteness
- Check Cartesian and polar coordinate consistency
