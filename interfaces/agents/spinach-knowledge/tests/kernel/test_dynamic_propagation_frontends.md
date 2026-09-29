# tests/kernel/test_dynamic_propagation_frontends.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_propagation_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_propagation_frontends.m)

## Purpose

Regression test for the dynamic propagation front-end kernels of Spinach. The test exercises `propagator()`, `step()`, `evolution()`, `krylov()`, and `reduce()` against direct finite-dimensional propagation references on tiny spin systems.

## Behaviour

The function announces the test target, initialises a regression test result via `new_test_result()` for the target `kernel/dynamic_propagation_frontends`, and then runs four sub-tests:

1. **`local_test_propagator_step`** — builds a one-spin spherical-tensor Liouville-space system (`sys.magnet = 14.1`, isotope `1H`, zero scalar Zeeman interaction, formalism `sphten-liouv`, approximation `none`, `assume(...,'nmr')`), forces Taylor propagation by setting `spin_system.tols.small_matrix = 2`, and uses `L = operator(spin_system,'Lz','1H')` with `rho = state(spin_system,'Lx','1H') + 0.25*state(spin_system,'Ly','1H')` and `dt = 2.5e-4`. It checks:
   - `propagator()` against `expm(full(-1i*L*dt))` (tolerances `1e-10`).
   - Numeric `step()` against the reference propagator action `P_ref*rho` (tolerances `1e-12`).
   - The zero-time shortcut: `step(spin_system,L,rho,0)` must return `rho` unchanged (tolerances `1e-14`).
   - Two-point (`{L,L}`) and three-point (`{L,L,L}`) quadrature paths reducing to the constant-generator result (tolerances `1e-12`).
   - The state-dependent generator dispatch `step(spin_system,{@(time,state)L,0.0,'PWCL'},rho,dt)` routed through `iserstep()` (tolerances `1e-12`).
   - The Hilbert-space commutator-series branch: a one-spin `zeeman-hilb` system with `spin_system.tols.small_matrix = 1`, `S = pauli(2)`, `H = 2*pi*123*S.z`, `rho_h = S.x + 0.25*S.y`, compared against `P*rho_h*P'` with `P = expm(-1i*H*dt)` (tolerances `1e-12`).

2. **`local_test_evolution`** — uses the Liouville system with trajectory-level reduction disabled (`sys.disable` augmented with `trajlevel`), `nsteps = 3`, and `P = propagator(spin_system,L,dt)`, building the reference trajectory `traj_ref = [rho P*rho P^2*rho P^3*rho]`. It checks the `final`, `trajectory`, `observable` (single coil `coil_x`), `multichannel` (coils `[coil_x coil_y]`), and `refocus` output modes of `evolution()` against explicit propagator products (tolerances `1e-12`). The refocus check propagates a stack `[rho 2*rho 3*rho]` for 2 steps so the nth stack member is propagated for n-1 steps. It also checks the `total` integral mode using the scalar damping generator `L_decay = -1i*speye(size(L,1))`, comparing against `real(coil_x'*rho)` (tolerances `1e-12`).

3. **`local_test_krylov`** — uses the same Liouville system, `dt`, `nsteps`, propagator, and reference trajectory to check the `final`, `trajectory`, `observable`, `multichannel`, and `refocus` output modes of `krylov()` against the same references (tolerances `1e-12`). The refocus call passes an empty `dt` (`krylov(spin_system,L,[],rho_stack,dt,[],'refocus')`).

4. **`local_test_reduce`** — calls `reduce(spin_system,L,rho)` with `rho = state(spin_system,'Lx','1H')` and checks that every returned projector has orthonormal columns (`P'*P` equals identity, tolerances `1e-14`), that the sum of projected components `sum_n P_n*(P_n'*rho)` reconstructs the source state (tolerances `1e-14`), and that with `trajlevel` disabled the function returns a scalar cell containing the identity.

All comparisons are registered through `test_close()` (and `test_true()` for the disable branch), accumulating messages into the returned regression test result.

## Inputs and outputs

```matlab
result = test_dynamic_propagation_frontends()
```

- **Output:** `result` — regression test result structure with explanatory messages for each checked branch.
- **Input:** none.

## References

- Spinach dynamic propagation kernels: `propagator()`, `step()`, `evolution()`, `krylov()`, `reduce()`.
- [Spinach GitHub repository](https://github.com/IlyaKuprov/Spinach)
