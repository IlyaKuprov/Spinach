# tests/kernel/test_states_composite_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_states_composite_suite.m`
- Signature: `result=test_states_composite_suite()`
- Total lines: 125

## Purpose

Tests composite state generators in kernel/states. Syntax: result=test_states_composite_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Composite state generators\n')`.
- Lines 19-22: State the state-generation target of the test; implemented by `result=new_test_result('kernel/states_composite_suite', 'Composite state generators', 'state helper functions must produce the textbook density operators.')`.
- Lines 24-25: Build a one-proton Hilbert-space spin system; implemented by `sys.magnet=0`.
- Lines 32-34: Check the Hilbert-space thermodynamic unit state; implemented by `result=test_close(result,'unit_state zeeman-hilb',unit_state(spin_system),speye(2),1e-15,1e-15, 'the Hilbert-space unit state is the identity density matrix')`.
- Lines 36-37: Check the Zeeman-Liouville unit vector normalisation; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 43-44: Build a two-proton Hilbert-space spin system; implemented by `clear sys inter bas`.
- Lines 54-55: Textbook two-spin wavefunctions in the Zeeman product basis; implemented by `alpha=[1;0]`.
- Lines 64-66: Check singlet and triplet projectors; implemented by `result=test_close(result,'singlet projector',singlet(spin_system,1,2),sing*sing',1e-14,1e-14, 'the singlet state is |alpha beta-beta alpha>< |/2')`.
- Lines 75-76: Build a four-proton Hilbert-space spin system; implemented by `clear sys inter bas`.
- Lines 88-89: Check the product of two singlet projectors; implemented by `singlet_pair=sing*sing'`.
- Lines 94-95: Build a three-proton Hilbert-space spin system for partner-state enumeration; implemented by `clear sys inter bas`.
- Lines 106-107: Enumerate two partner spins that may each be E or Lz while spin 2 is L+; implemented by `[A,descr]=partner_state(spin_system,{{'L+',2}},{{{'E','Lz'},[1 3]}})`.
- Lines 115-116: Check every enumerated partner state against a direct state() call; implemented by `full_spin_list={1,2,3}`.

### Control flow inferred from the code

- Line 117: `for` loop over `n=1:numel(expected_descr)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/states_composite_suite', 'Composite state generators', 'state helper functions must produce the textbook density operators.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 39: computes `unit` using `unit=speye(2); unit=unit(:); unit=unit/norm(unit,2)`.
- Lines 48: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=0`.
- Lines 49: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 55: computes `alpha` using `alpha=[1;0]`.
- Lines 56: computes `beta` using `beta=[0;1]`.
- Lines 57: computes `aa` using `aa=kron(alpha,alpha)`.
- Lines 58: computes `ab` using `ab=kron(alpha,beta)`.
- Lines 59: computes `ba` using `ba=kron(beta,alpha)`.
- Lines 60: computes `bb` using `bb=kron(beta,beta)`.
- Lines 61: computes `sing` using `sing=(ab-ba)/sqrt(2)`.
- Lines 62: computes `trip_zero` using `trip_zero=(ab+ba)/sqrt(2)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks unit-state normalisation, two-spin singlet/triplet
- projectors, four-spin product states, and partner-state enumeration.

## Implementation structure

- Tests composite state generators in kernel/states. Syntax:
- result=test_states_composite_suite()
- result -regression test result with explanatory messages
- The test checks unit-state normalisation, two-spin singlet/triplet
- projectors, four-spin product states, and partner-state enumeration.
- Announce the test target
- State the state-generation target of the test
- Build a one-proton Hilbert-space spin system
- Check the Hilbert-space thermodynamic unit state
- Check the Zeeman-Liouville unit vector normalisation
- Build a two-proton Hilbert-space spin system
- Textbook two-spin wavefunctions in the Zeeman product basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `test_close()`, `unit_state()`, `speye()`, `unit()`, `singlet()`, `triplet()`, `four_spin_states()`, `partner_state()`, `test_true()`, `isequal()`, `descr()`, `expected_descr()`, `state()`, `int2str()`.
