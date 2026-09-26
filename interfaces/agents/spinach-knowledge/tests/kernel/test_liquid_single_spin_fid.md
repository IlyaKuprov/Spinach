# tests/kernel/test_liquid_single_spin_fid.m

- Signature: `result=test_liquid_single_spin_fid()`

## Purpose

Tests a one-spin liquid-state free induction decay. Syntax: result=test_liquid_single_spin_fid()

## Physical / mathematical content

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
