# tests/kernel/test_state_constructor_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_state_constructor_suite.m`
- Signature: `result=test_state_constructor_suite()`
- Total lines: 120

## Purpose

Tests state-constructor helper functions. Syntax: result=test_state_constructor_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: State-constructor functions\n')`.
- Lines 20-23: State the state-constructor target of the test; implemented by `result=new_test_result('kernel/state_constructor_suite', 'State-constructor functions', 'state constructors must produce correctly normalised physical density objects.')`.
- Lines 25-26: One-spin Hilbert and Liouville unit states have known forms; implemented by `sys.magnet=0`.
- Lines 50-51: Two-spin singlet and triplet constructors must form four orthogonal projectors summing to identity; implemented by `sys2.magnet=0`.
- Lines 69-70: partner_state must enumerate all requested partner-state combinations; implemented by `sys3.magnet=0`.
- Lines 82-83: Four-spin singlet-singlet state must match the explicit product of two two-spin singlets; implemented by `sys4.magnet=0`.
- Lines 93-94: Spin-1 pair states must resolve the nine-dimensional Hilbert identity; implemented by `sysd.magnet=0`.
- Lines 107-108: Zero-field triplet projection must preserve density-matrix trace and Hermiticity; implemented by `syse.magnet=0`.

### Control flow inferred from the code

- Line 77: `for` loop over `n=1:numel(A)`.
- Line 100: `for` loop over `n=1:numel(Td), projector_sum=projector_sum+Td{n}; end`.
- Line 101: `for` loop over `n=1:numel(Qd), projector_sum=projector_sum+Qd{n}; end`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/state_constructor_suite', 'State-constructor functions', 'state constructors must produce correctly normalised physical density objects.')`.
- Lines 26: computes `sys.magnet` using `sys.magnet=0`.
- Lines 27: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 29: computes `inter.temperature` using `inter.temperature=300`.
- Lines 30: computes `bas.formalism` using `bas.formalism='zeeman-hilb'; bas.approximation='none'`.
- Lines 31: computes `spin_h` using `spin_h=test_spin_system(sys,inter,bas)`.
- Lines 37: computes `spin_l` using `spin_l=test_spin_system(sys,inter,bas)`.
- Lines 38: computes `unit_zeeman` using `unit_zeeman=speye(2); unit_zeeman=unit_zeeman(:)/sqrt(2)`.
- Lines 42: computes `spin_s` using `spin_s=test_spin_system(sys,inter,bas)`.
- Lines 43: computes `unit_s` using `unit_s=unit_state(spin_s)`.
- Lines 51: computes `sys2.magnet` using `sys2.magnet=0`.
- Lines 52: computes `sys2.isotopes` using `sys2.isotopes={'1H','1H'}`.
- Lines 53: computes `inter2.zeeman.scalar` using `inter2.zeeman.scalar={0,0}`.
- Lines 54: computes `bas2.formalism` using `bas2.formalism='zeeman-hilb'; bas2.approximation='none'`.
- Lines 55: computes `spin2` using `spin2=test_spin_system(sys2,inter2,bas2)`.
- Lines 56: computes `S` using `S=singlet(spin2,1,2)`.
- Lines 57: computes `[TU,T0,TD]` using `[TU,T0,TD]=triplet(spin2,1,2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks unit, thermal, singlet/triplet, partner, deuteron-pair,
- four-spin, and zero-field triplet state constructors using projector and
- normalisation identities.

## Implementation structure

- Tests state-constructor helper functions. Syntax:
- result=test_state_constructor_suite()
- result -regression test result with explanatory messages
- The test checks unit, thermal, singlet/triplet, partner, deuteron-pair,
- four-spin, and zero-field triplet state constructors using projector and
- normalisation identities.
- Announce the test target
- State the state-constructor target of the test
- One-spin Hilbert and Liouville unit states have known forms
- Two-spin singlet and triplet constructors must form four orthogonal projectors summing to identity
- partner_state must enumerate all requested partner-state combinations
- Four-spin singlet-singlet state must match the explicit product of two two-spin singlets

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `test_close()`, `unit_state()`, `speye()`, `equilibrium()`, `unit_zeeman()`, `stateinfo()`, `test_true()`, `singlet()`, `triplet()`, `partner_state()`, `int2str()`, `state()`, `four_spin_states()`, `deut_pair()`.
