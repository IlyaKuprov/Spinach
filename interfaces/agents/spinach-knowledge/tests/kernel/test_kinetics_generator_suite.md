# tests/kernel/test_kinetics_generator_suite.m

- Signature: `result=test_kinetics_generator_suite()`

## Purpose

Tests kinetics and flow generator helpers. Syntax: result=test_kinetics_generator_suite()

## Physical / mathematical content

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
