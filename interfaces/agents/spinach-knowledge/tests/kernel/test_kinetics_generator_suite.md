# tests/kernel/test_kinetics_generator_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_kinetics_generator_suite.m`
- Signature: `result=test_kinetics_generator_suite()`
- Total lines: 76

## Purpose

Tests kinetics and flow generator helpers. Syntax: result=test_kinetics_generator_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Kinetics and flow generator functions\n')`.
- Lines 20-23: State the kinetics target of the test; implemented by `result=new_test_result('kernel/kinetics_generator_suite', 'Kinetics and flow generator functions', 'kinetic generators must conserve matter and equilibrate closed system…`.
- Lines 25-26: A two-state reversible Markov generator has a closed equilibrium ratio; implemented by `K=[-2 1;2 -1]`.
- Lines 34-35: Build a two-site exchange Spinach system used by react_gen and kinetics; implemented by `sys.magnet=0`.
- Lines 44-45: react_gen must build a conservative drain/fill mapping for a specified reaction; implemented by `reaction.reactants=1`.
- Lines 54-55: Full kinetics generator for symmetric exchange must conserve population; implemented by `Kspin=kinetics(spin_system)`.
- Lines 59-60: A minimal two-cell diffusion mesh must produce a conservative symmetric generator; implemented by `mesh.vor.ncells=2`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/kinetics_generator_suite', 'Kinetics and flow generator functions', 'kinetic generators must conserve matter and equilibrate closed system…`.
- Lines 26: computes `K` using `K=[-2 1;2 -1]`.
- Lines 27: computes `c0` using `c0=[3;0]`.
- Lines 28: computes `ceq` using `ceq=equilibrate(K,c0)`.
- Lines 30: computes `'at equilibrium k_21*c_1` using `'at equilibrium k_21*c_1=k_12*c_2 while total concentration is conserved')`.
- Lines 35: computes `sys.magnet` using `sys.magnet=0`.
- Lines 36: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 37: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0}`.
- Lines 38: computes `inter.chem.parts` using `inter.chem.parts={1,2}`.
- Lines 39: computes `inter.chem.rates` using `inter.chem.rates=[-1 1;1 -1]`.
- Lines 40: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.
- Lines 41: computes `bas.formalism` using `bas.formalism='sphten-liouv'; bas.approximation='none'`.
- Lines 42: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 45: computes `reaction.reactants` using `reaction.reactants=1`.
- Lines 46: computes `reaction.products` using `reaction.products=2`.
- Lines 47: computes `reaction.matching` using `reaction.matching=[1 2]`.
- Lines 48: computes `G` using `G=react_gen(spin_system,reaction)`.
- Lines 55: computes `Kspin` using `Kspin=kinetics(spin_system)`.

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
