# tests/kernel/test_slr_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_slr_pulse.m`
- Signature: `result=test_slr_pulse()`
- Total lines: 191

## Purpose

Tests Shinnar-Le Roux selective excitation pulse design. Syntax: result=test_slr_pulse()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `throws_with()`, `ck_profile()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `slr_pulse()`, `test_true()`, `isequal()`, `all()`, `test_close()`, `durs()`, `small_x()`, `small_y()`, `small_durs()`, `num2str()`, `ck_profile()`, `log10()`, `transverse()`, `longitudinal()`, `small_trans()`.
