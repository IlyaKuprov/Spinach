# tests/kernel/test_scalar_coupling_hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_scalar_coupling_hamiltonian.m`
- Signature: `result=test_scalar_coupling_hamiltonian()`
- Total lines: 46

## Purpose

Tests the two-spin scalar-coupling Hamiltonian. Syntax: result=test_scalar_coupling_hamiltonian()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Scalar coupling Hamiltonian\n')`.
- Lines 19-22: State the Hamiltonian target of the test; implemented by `result=new_test_result('kernel/scalar_coupling_hamiltonian', 'Scalar coupling Hamiltonian', 'an isotropic J coupling must produce 2*pi*J I dot S.')`.
- Lines 24-25: Build a two-proton Hilbert-space spin system with a 10 Hz J coupling; implemented by `sys.magnet=0`.
- Lines 34-35: Build Spinach and textbook Hamiltonians; implemented by `H_obs=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 41-43: Check the scalar-coupling Hamiltonian; implemented by `result=test_close(result,'isotropic J Hamiltonian',H_obs,H_ref,1e-9,1e-12, 'scalar coupling is rotationally invariant I dot S in rad/s units')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/scalar_coupling_hamiltonian', 'Scalar coupling Hamiltonian', 'an isotropic J coupling must produce 2*pi*J I dot S.')`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0`.
- Lines 26: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.
- Lines 28: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 29: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 30: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=test_spin_system(sys,inter,bas)`.
- Lines 35: computes `H_obs` using `H_obs=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 36: computes `IxSx` using `IxSx=operator(spin_system,{'Lx','Lx'},{1,2})`.
- Lines 37: computes `IySy` using `IySy=operator(spin_system,{'Ly','Ly'},{1,2})`.
- Lines 38: computes `IzSz` using `IzSz=operator(spin_system,{'Lz','Lz'},{1,2})`.
- Lines 39: computes `H_ref` using `H_ref=2*pi*10*(IxSx+IySy+IzSz)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that an isotropic scalar coupling J produces the textbook
- Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).

## Implementation structure

- Tests the two-spin scalar-coupling Hamiltonian. Syntax:
- result=test_scalar_coupling_hamiltonian()
- result -regression test result with explanatory messages
- The test checks that an isotropic scalar coupling J produces the textbook
- Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).
- Announce the test target
- State the Hamiltonian target of the test
- Build a two-proton Hilbert-space spin system with a 10 Hz J coupling
- Build Spinach and textbook Hamiltonians
- Check the scalar-coupling Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `hamiltonian()`, `assume()`, `operator()`, `test_close()`.
