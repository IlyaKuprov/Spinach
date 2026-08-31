# tests/kernel/test_zeeman_hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_zeeman_hamiltonian.m`
- Signature: `result=test_zeeman_hamiltonian()`
- Total lines: 42

## Purpose

Tests the one-spin Zeeman Hamiltonian. Syntax: result=test_zeeman_hamiltonian()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Zeeman Hamiltonian sign and units\n')`.
- Lines 19-22: State the Hamiltonian target of the test; implemented by `result=new_test_result('kernel/zeeman_hamiltonian', 'Zeeman Hamiltonian sign and units', 'a scalar chemical shift must enter the NMR Hamiltonian with Spinach sign and ra…`.
- Lines 24-25: Build a one-proton Hilbert-space spin system with a 1 ppm shift; implemented by `sys.magnet=14.1`.
- Lines 32-33: Build Spinach and reference Hamiltonians; implemented by `H_obs=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 37-39: Check the physical frequency and sign convention; implemented by `result=test_close(result,'one-spin Zeeman Hamiltonian',H_obs,H_ref,1e-6,1e-12, 'positive ppm gives -omega*Lz in the Spinach NMR rotating-frame convention')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/zeeman_hamiltonian', 'Zeeman Hamiltonian sign and units', 'a scalar chemical shift must enter the NMR Hamiltonian with Spinach sign and ra…`.
- Lines 25: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 33: computes `H_obs` using `H_obs=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 34: computes `nu` using `nu=ppm2hz(1,sys.magnet,'1H')`.
- Lines 35: computes `H_ref` using `H_ref=-2*pi*nu*operator(spin_system,'Lz',1)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks Spinach's NMR convention for a positive chemical shift:
- the rotating-frame Hamiltonian contribution is -2*pi*nu*Lz.

## Implementation structure

- Tests the one-spin Zeeman Hamiltonian. Syntax:
- result=test_zeeman_hamiltonian()
- result -regression test result with explanatory messages
- The test checks Spinach's NMR convention for a positive chemical shift:
- the rotating-frame Hamiltonian contribution is -2*pi*nu*Lz.
- Announce the test target
- State the Hamiltonian target of the test
- Build a one-proton Hilbert-space spin system with a 1 ppm shift
- Build Spinach and reference Hamiltonians
- Check the physical frequency and sign convention

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `hamiltonian()`, `assume()`, `ppm2hz()`, `operator()`, `test_close()`.
