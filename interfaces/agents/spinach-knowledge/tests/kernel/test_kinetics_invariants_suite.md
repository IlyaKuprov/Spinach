# tests/kernel/test_kinetics_invariants_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_kinetics_invariants_suite.m`
- Signature: `result=test_kinetics_invariants_suite()`
- Total lines: 79

## Purpose

Tests deterministic chemical kinetics helpers. Syntax: result=test_kinetics_invariants_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks closed-form steady states, independent reaction blocks,
- reaction-generator state routing, and conservation in a tiny exchange
- kinetics superoperator.

## Implementation structure

- Tests deterministic chemical kinetics helpers. Syntax:
- result=test_kinetics_invariants_suite()
- result -regression test result with explanatory messages
- The test checks closed-form steady states, independent reaction blocks,
- reaction-generator state routing, and conservation in a tiny exchange
- kinetics superoperator.
- Announce the test target
- State the kinetics target of the test
- Check a two-site steady state from detailed balance
- Check recursive treatment of independent reaction blocks
- Check the zero-concentration shortcut
- Build a two-site spherical-tensor spin system for exchange and reactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `equilibrate()`, `blkdiag()`, `test_spin_system()`, `kinetics()`, `react_gen()`, `state()`.
