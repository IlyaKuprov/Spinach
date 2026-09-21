# tests/kernel/test_dynamic_state_projection_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_state_projection_suite.m`
- Signature: `result=test_dynamic_state_projection_suite()`
- Total lines: 110

## Purpose

Tests dynamic state projection helper paths. Syntax: result=test_dynamic_state_projection_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks deuteron-pair coherences, dephased population stationarity,
- captured stateinfo() output, and isotropic zero-field triplet projection.

## Implementation structure

- Tests dynamic state projection helper paths. Syntax:
- result=test_dynamic_state_projection_suite()
- result -regression test result with explanatory messages
- The test checks deuteron-pair coherences, dephased population stationarity,
- captured stateinfo() output, and isotropic zero-field triplet projection.
- Announce the test target
- State the state-projection target of the test
- Build a two-deuteron Hilbert-space system
- Request all deuteron-pair coherences
- Check adjoint pairing of triplet coherences
- Check adjoint pairing of quintet coherences
- Request dephased deuteron-pair populations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `deut_pair()`, `test_close()`, `operator()`, `int2str()`, `evalc()`, `stateinfo()`, `state()`, `test_true()`, `contains()`, `report()`, `zftrip()`, `speye()`.
