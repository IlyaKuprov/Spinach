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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Chemical kinetics invariants\n')`.
- Lines 20-23: State the kinetics target of the test; implemented by `result=new_test_result('kernel/kinetics_invariants_suite', 'Chemical kinetics invariants', 'kinetic generators must conserve matter and route spin order between declared…`.
- Lines 25-26: Check a two-site steady state from detailed balance; implemented by `kf=2`.
- Lines 35-36: Check recursive treatment of independent reaction blocks; implemented by `K1=[-1 4;1 -4]`.
- Lines 44-46: Check the zero-concentration shortcut; implemented by `result=test_close(result,'equilibrate zero concentration',equilibrate(K,zeros(4,1)),zeros(4,1),1e-15,1e-15, 'a zero initial concentration vector remains zero for linear…`.
- Lines 48-49: Build a two-site spherical-tensor spin system for exchange and reactions; implemented by `sys.magnet=14.1`.
- Lines 59-60: Closed exchange must conserve every column sum of the kinetic generator; implemented by `K=kinetics(spin_system)`.
- Lines 65-66: A one-way reaction generator must drain source spin order and fill matched product spin order; implemented by `reaction.reactants=1`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/kinetics_invariants_suite', 'Chemical kinetics invariants', 'kinetic generators must conserve matter and route spin order between declared…`.
- Lines 26: computes `kf` using `kf=2`.
- Lines 27: computes `kr` using `kr=5`.
- Lines 28: computes `K` using `K=[-kf kr; kf -kr]`.
- Lines 29: computes `c0` using `c0=[2;1]`.
- Lines 30: computes `ctot` using `ctot=sum(c0)`.
- Lines 31: computes `c_ref` using `c_ref=ctot*[kr; kf]/(kf+kr)`.
- Lines 36: computes `K1` using `K1=[-1 4;1 -4]`.
- Lines 37: computes `K2` using `K2=[-3 2;3 -2]`.
- Lines 49: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 50: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 51: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.
- Lines 52: computes `inter.chem.parts` using `inter.chem.parts={1,2}`.
- Lines 53: computes `inter.chem.rates` using `inter.chem.rates=[-3 3;3 -3]`.
- Lines 54: computes `inter.chem.concs` using `inter.chem.concs=[1 1]`.
- Lines 55: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 56: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 57: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.

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
