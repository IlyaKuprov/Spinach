# tests/kernel/test_kinetics_generator_suite.m

- Signature: `result=test_kinetics_generator_suite()`

## Purpose

Tests kinetics and flow generator helpers for equilibrium, population conservation, and a minimal diffusion case.

## Assertions

- For `K=[-2 1;2 -1]` and initial population `[3;0]`, `equilibrate` returns `[1;2]`, conserving total population and satisfying the stationary balance of the two-state generator. With zero initial population, it returns `[0;0]` exactly.
- Constructs a two-site exchange Spinach system with isotopes `{'1H','1H'}`, zero Zeeman scalars, chemical parts `{1,2}`, rates `[-1 1;1 -1]`, concentrations `[1 1]`, and a `sphten-liouv` basis with approximation `none`. For the reaction with reactant `1`, product `2`, and matching `[1 2]`, `react_gen` returns one generator matrix whose columns sum to zero within `1e-14` absolute and relative tolerances.
- The full `kinetics` generator for that exchange system has zero column sums within `1e-14` absolute and relative tolerances.
- Builds a two-cell diffusion mesh with unit cell weights, vertices `[0 0;0 1]`, both cells using `[1 2]`, active indices `[1;2]`, triangle `[1 2 3]`, coordinates `x=[0;1;0]` and `y=[0;0;1]`, and zero `u` and `v` fields. With `diff=0.5`, `flow_gen` returns `[-0.5 0.5;0.5 -0.5]` and has column sums `[0 0]`, each within `1e-14` absolute and relative tolerances.

## Output

- `result` — regression test result with explanatory messages.