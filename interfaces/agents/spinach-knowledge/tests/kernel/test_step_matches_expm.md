# tests/kernel/test_step_matches_expm.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_matches_expm.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_matches_expm.m)

## Purpose

Regression test that verifies Hilbert-space propagation performed by `step()` against direct matrix exponentiation with `expm()`. The test checks the Spinach sign convention for density-matrix evolution:

`rho(t) = exp(-iHt) rho(0) exp(+iHt)`

## Behaviour

- Announces the test target with `fprintf('TESTING: Hilbert propagation against expm\n')`.
- Initialises a regression test result via `new_test_result('kernel/step_matches_expm', 'Hilbert propagation against expm', 'step() must reproduce unitary density-matrix propagation.')`.
- Builds a one-proton Hilbert-space spin system with `sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.formalism='zeeman-hilb'`, and `bas.approximation='none'`, using `test_spin_system(sys,inter,bas)`.
- Defines the Hamiltonian `H = 2*pi*123*S.z` and the initial density matrix `rho = S.x + 0.25*S.y` from Pauli matrices `S=pauli(2)`, with a time step `dt = 2.5e-3`.
- Constructs the independent exact propagator `P = expm(-1i*H*dt)` and computes the reference evolved state `rho_ref = P*rho*P'`.
- Runs `rho_obs = step(spin_system,H,rho,dt)` and compares the two with `test_close(result,'step versus expm',rho_obs,rho_ref,1e-13,1e-13,'finite Hilbert-space propagation is exactly unitary')`.

## Inputs and outputs

```matlab
result = test_step_matches_expm()
```

- **Output:** `result` — regression test result with explanatory messages.
- **Input:** none.

## References

- Source file: [tests/kernel/test_step_matches_expm.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_matches_expm.m) in the Spinach repository.
