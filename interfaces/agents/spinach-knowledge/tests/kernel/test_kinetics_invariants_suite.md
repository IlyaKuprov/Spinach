# tests/kernel/test_kinetics_invariants_suite.m

- Signature: `result=test_kinetics_invariants_suite()`

## Purpose

Tests deterministic chemical kinetics helpers. Syntax: result=test_kinetics_invariants_suite()

## Physical / mathematical content

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
