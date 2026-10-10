# tests/kernel/test_step_zero_time.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_zero_time.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_zero_time.m)

## Purpose

Regression test for zero-duration propagation. The test verifies the identity limit of the propagator: a zero time step must leave the density matrix exactly unchanged.

## Behaviour

- Announces the test target with `fprintf('TESTING: Zero-duration propagation identity\n')`.
- Registers a new test result via `new_test_result` with test name `'kernel/step_zero_time'`, description `'Zero-duration propagation identity'`, and the statement `'a propagator over zero time is the identity map.'`.
- Builds a one-proton Hilbert-space spin system with `sys.magnet=0`, `sys.isotopes={'1H'}`, `inter.zeeman.scalar={0}`, `bas.formalism='zeeman-hilb'`, and `bas.approximation='none'`, using `test_spin_system(sys,inter,bas)`.
- Constructs an arbitrary Hermitian density matrix `rho=S.x+2*S.z` and a Hamiltonian `H=3*S.x+5*S.z` from `pauli(2)`.
- Propagates the density matrix for zero time with `step(spin_system,H,rho,0)`.
- Checks the identity limit with `test_close(result,'rho(t=0)=rho(0)',rho_obs,rho,1e-15,1e-15,'zero-duration evolution cannot change the state')`, using absolute and relative tolerances of `1e-15`.

## Inputs and outputs

```matlab
result=test_step_zero_time()
```

- **Output:** `result` — regression test result with explanatory messages.
- **Input:** none.

## References

- [tests/kernel/test_step_zero_time.m — Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_step_zero_time.m)
