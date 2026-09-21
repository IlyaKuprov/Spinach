# tests/kernel/test_liquid_single_spin_fid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_liquid_single_spin_fid.m`
- Signature: `result=test_liquid_single_spin_fid()`
- Total lines: 52

## Purpose

Tests a one-spin liquid-state free induction decay. Syntax: result=test_liquid_single_spin_fid()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test simulates a zero-offset one-spin FID. With no Hamiltonian and no
- relaxation, transverse magnetisation is constant in time.

## Implementation structure

- Tests a one-spin liquid-state free induction decay. Syntax:
- result=test_liquid_single_spin_fid()
- result -regression test result with explanatory messages
- The test simulates a zero-offset one-spin FID. With no Hamiltonian and no
- relaxation, transverse magnetisation is constant in time.
- Announce the test target
- State the NMR target of the test
- Build a one-spin Liouville-space system
- Set up a zero-offset acquisition
- Simulate the FID
- The signal is constant and equal to its first point

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `state()`, `liquid()`, `test_close()`, `fid()`.
