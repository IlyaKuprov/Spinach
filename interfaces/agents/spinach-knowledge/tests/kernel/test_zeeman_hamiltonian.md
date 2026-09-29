# tests/kernel/test_zeeman_hamiltonian.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_zeeman_hamiltonian.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_zeeman_hamiltonian.m)

## Purpose

Regression test for the one-spin Zeeman Hamiltonian. The test verifies Spinach's NMR convention for a positive chemical shift: the rotating-frame Hamiltonian contribution is `-2*pi*nu*Lz`, i.e. a positive ppm shift enters the NMR Hamiltonian with the Spinach sign convention and rad/s units.

## Behaviour

- Announces the test target with `fprintf('TESTING: Zeeman Hamiltonian sign and units\n')`.
- Initialises a regression test result via `new_test_result('kernel/zeeman_hamiltonian', 'Zeeman Hamiltonian sign and units', 'a scalar chemical shift must enter the NMR Hamiltonian with Spinach sign and rad/s units.')`.
- Builds a one-proton Hilbert-space spin system with a 1 ppm scalar shift:
  - `sys.magnet = 14.1`
  - `sys.isotopes = {'1H'}`
  - `inter.zeeman.scalar = {1}`
  - `bas.formalism = 'zeeman-hilb'`
  - `bas.approximation = 'none'`
  - the system is assembled with `test_spin_system(sys, inter, bas)`.
- Constructs the Spinach Hamiltonian as `H_obs = hamiltonian(assume(spin_system, 'nmr'))`.
- Computes the reference frequency `nu = ppm2hz(1, sys.magnet, '1H')` and the reference Hamiltonian `H_ref = -2*pi*nu*operator(spin_system, 'Lz', 1)`.
- Compares `H_obs` and `H_ref` with `test_close(result, 'one-spin Zeeman Hamiltonian', H_obs, H_ref, 1e-6, 1e-12, 'positive ppm gives -omega*Lz in the Spinach NMR rotating-frame convention')`, using absolute tolerance `1e-6` and relative tolerance `1e-12`.

## Inputs and outputs

- **Inputs:** none. The function takes no arguments.
- **Outputs:**
  - `result` — regression test result structure with explanatory messages, as produced by `new_test_result` and updated by `test_close`.

## References

- [Spinach — tests/kernel/test_zeeman_hamiltonian.m (GitHub)](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_zeeman_hamiltonian.m)
