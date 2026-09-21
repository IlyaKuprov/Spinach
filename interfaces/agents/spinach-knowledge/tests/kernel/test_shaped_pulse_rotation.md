# tests/kernel/test_shaped_pulse_rotation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_shaped_pulse_rotation.m`
- Signature: `result=test_shaped_pulse_rotation()`
- Total lines: 49

## Purpose

Tests a one-slice Cartesian shaped pulse. Syntax: result=test_shaped_pulse_rotation()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test applies one rectangular X pulse slice with amplitude 1 rad/s and
- duration pi seconds; the net flip angle is pi, so Lz must invert.

## Implementation structure

- Tests a one-slice Cartesian shaped pulse. Syntax:
- result=test_shaped_pulse_rotation()
- result -regression test result with explanatory messages
- The test applies one rectangular X pulse slice with amplitude 1 rad/s and
- duration pi seconds; the net flip angle is pi, so Lz must invert.
- Announce the test target
- State the pulse target of the test
- Build a one-proton Hilbert-space spin system
- Define drift, controls, and one pi pulse slice
- Apply the shaped pulse
- Check the physical rotation result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `operator()`, `state()`, `shaped_pulse_xy()`, `test_close()`.
