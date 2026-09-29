# tests/kernel/test_dynamic_equilibrium_frontends.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_equilibrium_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_equilibrium_frontends.m)

## Purpose

Regression test for the dynamic equilibrium front-end kernels: `thermalize()`, `steady()`, and `residual()`. The test exercises these functions on compact Liouville-space systems with explicit fixed-point references, verifying that thermalisation and steady-state helpers produce explicit fixed points.

## Behaviour

- Announces the test target with `TESTING: Dynamic equilibrium front ends` and registers a test result under `kernel/dynamic_equilibrium_frontends` via `new_test_result()`.
- Runs three local subtests:
  - `local_test_thermalize` — checks IME and DiBari thermalisation branches.
  - `local_test_steady` — checks Newton and repeated-squaring steady-state solvers.
  - `local_test_residual` — checks weak residual-order tensor reduction.
- **Thermalise subtest:** builds a one-spin spherical-tensor Liouville-space system (`local_liouville_system`), forms a sparse relaxation superoperator `R = -diag([0; ones(dim-1,1)])` and an equilibrium state `rho_eq = unit + 0.1*state(spin_system,'Lz','1H')`.
  - IME mode: `thermalize(spin_system,R,[],[],rho_eq,'IME')` must make the requested equilibrium state stationary, i.e. `R_ime*rho_eq` matches a zero vector to tolerances `1e-14` (absolute and relative).
  - DiBari mode: `thermalize(spin_system,R,H_left,temperature,[],'dibari')` with `H_left = operator(spin_system,'Lz','1H','left')` and `temperature = 300.0` must equal the reference product `R*propagator(spin_system,H_left,1i*beta)`, where `beta = spin_system.tols.hbar/(spin_system.tols.kbol*temperature)`, to tolerances `1e-14`.
- **Steady subtest:** builds the same one-spin system and constructs a contractive affine propagator `P` with a known fixed point `rho_ss = [1.0; 0.20; -0.10; 0.05]` using `contract = diag([0.25 0.50 0.75])`:
  - Verifies `P*rho_ss` equals `rho_ss` to tolerances `1e-14`.
  - Newton mode: `steady(spin_system,P,[],'newton')` must recover the fixed point to tolerances `1e-12`.
  - Squaring mode: `steady(spin_system,P,[],'squaring')` must converge to the same fixed point to tolerances `1e-10`.
- **Residual subtest:** builds a heteronuclear `1H`/`13C` spin system (`sys.magnet = 5.9`) with Zeeman scalars `{5.0, 65.0}`, scalar coupling `140.0` between spins 1 and 2, coordinates `[0.0 0.0 0.0]` and `[0.6 0.7 0.8]`, an order matrix `diag([1e-3 2e-3 -3e-3])`, formalism `sphten-liouv`, and approximation `none`. It stores the coupling tensor before and after `residual(spin_system)` and checks:
  - `trace(J_after)` equals `trace(J_before)` to tolerances `1e-10`/`1e-12` (isotropic part preserved).
  - `J_after(1,1)` equals `J_after(2,2)` to tolerances `1e-12` (axially symmetric residual tensor).
  - Off-diagonal components of `J_after` are zero to tolerances `1e-12`.
- **Helper `local_liouville_system`:** builds a one-spin `1H` system with `sys.magnet = 14.1`, zero Zeeman scalar, formalism `sphten-liouv`, approximation `none`, via `test_spin_system()`, then applies `assume(spin_system,'nmr')`.

## Inputs and outputs

**Syntax:**

```matlab
result = test_dynamic_equilibrium_frontends()
```

**Outputs:**

- `result` — regression test result structure with explanatory messages, accumulated through `test_close()` checks.

## References
