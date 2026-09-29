# tests/kernel/test_dynamic_state_equilibrium_suite.m

## Purpose

Regression test suite for the thermal equilibrium state constructors in Spinach. The suite verifies that `equilibrium()` produces normalised Boltzmann states in the supported formalisms, checking Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian construction paths against direct Boltzmann references.

## Behaviour

- Announces the test target with `fprintf('TESTING: Thermal equilibrium constructors\n')` and initialises a regression test result via `new_test_result()` for `kernel/dynamic_state_equilibrium_suite`, with the specification that `equilibrium()` must produce normalised Boltzmann states in supported formalisms.
- Builds a one-spin Hilbert-space system using `test_spin_system()` with `sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `inter.temperature=300`, `bas.formalism='zeeman-hilb'`, and `bas.approximation='none'`.
- Sets an explicit non-degenerate Hamiltonian `H=2*pi*diag([-5e11 3e11])` in angular frequency units and computes the reference Boltzmann density matrix `rho_ref=expm(-beta*H)` normalised by its trace, where `beta=spin_h.tols.hbar/(spin_h.tols.kbol*spin_h.rlx.temperature)`.
- Hilbert-space checks:
  - `equilibrium(spin_h,H)` must equal `rho_ref` within absolute and relative tolerances of `1e-13`, with the message that Hilbert-space thermal equilibrium is `exp(-beta H)/Tr[exp(-beta H)]`.
  - `trace(rho_h)` must equal 1 within `1e-13`, since a density matrix returned by `equilibrium()` must have unit trace.
  - `rho_h` must equal its Hermitian conjugate `rho_h'` within `1e-13`, since a Boltzmann density matrix of a Hermitian Hamiltonian is Hermitian.
- Builds the matching Zeeman-Liouville system with `bas.formalism='zeeman-liouv'` and forms the left-product Hamiltonian `H_left=kron(speye(2),H)`; `equilibrium(spin_l,H_left)` must match the vectorised reference `rho_ref(:)` within `1e-13`, since left-product Liouville equilibrium must match the vectorised Hilbert density matrix.
- Oriented-Hamiltonian branch:
  - Constructs a minimal anisotropic Hamiltonian cell `Q{1}` as a 3-by-3 cell array of sparse 2-by-2 matrices, with `Q{1}{2,2}=sparse(2*pi*diag([2e11 -2e11]))` and all other entries sparse 2-by-2 zero matrices.
  - Uses `euler_angles=[0 0 0]` and computes `H_oriented=H+orientation(Q,euler_angles)`.
  - Compares the four-argument call `equilibrium(spin_h,H,Q,euler_angles)` against `equilibrium(spin_h,H_oriented)` within `1e-13`, with the message that the oriented branch must thermalise `H+orientation(Q,euler_angles)`.
- All comparisons are performed through `test_close()`, which accumulates results and explanatory messages into the returned test result.

## Inputs and outputs

```matlab
result = test_dynamic_state_equilibrium_suite()
```

- **Outputs**:
  - `result` — regression test result with explanatory messages.
- **Inputs**: none.

## References

- Source: [tests/kernel/test_dynamic_state_equilibrium_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_state_equilibrium_suite.m)
