# tests/kernel/test_kinetics_generator_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_kinetics_generator_suite.m`
- Signature: `result=test_kinetics_generator_suite()`
- Total lines: 76

## Purpose

Tests kinetics and flow generator helpers. Syntax: result=test_kinetics_generator_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks linear equilibrium, chemical reaction generators, full
- chemical kinetics generators, and a minimal hydrodynamic diffusion
- generator against conservation and detailed-balance invariants.

## Implementation structure

- Tests kinetics and flow generator helpers. Syntax:
- result=test_kinetics_generator_suite()
- result -regression test result with explanatory messages
- The test checks linear equilibrium, chemical reaction generators, full
- chemical kinetics generators, and a minimal hydrodynamic diffusion
- generator against conservation and detailed-balance invariants.
- Announce the test target
- State the kinetics target of the test
- A two-state reversible Markov generator has a closed equilibrium ratio
- Build a two-site exchange Spinach system used by react_gen and kinetics
- react_gen must build a conservative drain/fill mapping for a specified reaction
- Full kinetics generator for symmetric exchange must conserve population

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `equilibrate()`, `test_close()`, `test_spin_system()`, `react_gen()`, `kinetics()`, `flow_gen()`.
