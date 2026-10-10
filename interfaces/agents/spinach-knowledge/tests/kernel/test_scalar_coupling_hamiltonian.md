# tests/kernel/test_scalar_coupling_hamiltonian.m

Source: [tests/kernel/test_scalar_coupling_hamiltonian.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_scalar_coupling_hamiltonian.m)

## Purpose

Regression test for the two-spin scalar-coupling Hamiltonian. The test verifies that an isotropic scalar coupling J produces the textbook Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).

## Behaviour

- Announces the test target with `fprintf('TESTING: Scalar coupling Hamiltonian\n')`.
- Initialises a regression test result via `new_test_result` for `kernel/scalar_coupling_hamiltonian`, describing the target as "an isotropic J coupling must produce 2*pi*J I dot S."
- Builds a two-proton Hilbert-space spin system with a 10 Hz J coupling:
  - `sys.magnet=0`
  - `sys.isotopes={'1H','1H'}`
  - `inter.zeeman.scalar={0,0}`
  - `inter.coupling.scalar{1,2}=10`
  - `inter.coupling.scalar{2,2}=0`
  - `bas.formalism='zeeman-hilb'`
  - `bas.approximation='none'`
  - `spin_system=test_spin_system(sys,inter,bas)`
- Computes the Spinach Hamiltonian as `H_obs=hamiltonian(assume(spin_system,'nmr'))`.
- Constructs the reference Hamiltonian from Cartesian spin operators:
  - `IxSx=operator(spin_system,{'Lx','Lx'},{1,2})`
  - `IySy=operator(spin_system,{'Ly','Ly'},{1,2})`
  - `IzSz=operator(spin_system,{'Lz','Lz'},{1,2})`
  - `H_ref=2*pi*10*(IxSx+IySy+IzSz)`
- Compares observed and reference Hamiltonians with `test_close(result,'isotropic J Hamiltonian',H_obs,H_ref,1e-9,1e-12,...)`, with the message "scalar coupling is rotationally invariant I dot S in rad/s units".

## Inputs and outputs

- Syntax: `result=test_scalar_coupling_hamiltonian()`
- The function takes no inputs.
- `result` — regression test result with explanatory messages.

## References

- [Spinach library](https://spindynamics.org/)
- [Spinach on GitHub](https://github.com/IlyaKuprov/Spinach)
