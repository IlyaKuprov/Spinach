# tests/kernel/test_kinetics_invariants_suite.m

- Signature: `result=test_kinetics_invariants_suite()`

## Purpose

Regression-tests deterministic chemical kinetics helpers: equilibrium concentrations, independent reaction blocks, reaction-generator state routing, and conservation in a two-site exchange kinetics superoperator.

## Physical / mathematical content

- For a two-site reaction with forward rate `kf=2` and reverse rate `kr=5`, checks that equilibrium satisfies `kf*c_1=kr*c_2` while conserving total concentration.
- Checks that independent reaction blocks equilibrate separately, each retaining its own total concentration, and that zero initial concentrations remain zero.
- In a closed two-site exchange model, checks that every column of the kinetic generator sums to zero.
- Checks that a one-way reaction generator removes matched spin order from the reactant species and inserts it into the product species, without draining spin order already on the product species.

## Numerical / algorithmic content

- Compares `equilibrate(K,c0)` against closed-form steady states for a two-site system and two independent reaction blocks, using absolute and relative tolerances of `1e-13`; the zero-concentration check uses `1e-15`.
- Constructs a two-site `1H` spin system at `14.1` with equal concentrations and exchange rates `[-3 3;3 -3]`, using the `sphten-liouv` formalism with no approximation. Checks column sums of `kinetics(spin_system)` and state routing by `react_gen(spin_system,reaction)` with tolerances of `1e-14`.

## Outputs

- `result` — regression test result with explanatory messages.