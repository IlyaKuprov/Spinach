# tests/kernel/test_dynamic_trajectory_frontends.m

**Source**: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_trajectory_frontends.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_trajectory_frontends.m)

## Purpose

Regression test for the dynamic trajectory-analysis front-end kernels. It exercises the plotting branches of `trajan()` and the scoring branches of `trajsimil()` on a compact two-spin spherical-tensor trajectory, checking that both helpers expose deterministic branch outputs.

## Behavior

- Announces the test target with `fprintf` and initializes a regression test result via `new_test_result` under the identifier `kernel/dynamic_trajectory_frontends`, described as "Dynamic trajectory analysis front ends" with the requirement that "trajectory plotting and similarity helpers must expose deterministic branch outputs".
- Forces figures to be invisible during plotting checks by setting the groot default figure visibility to `'off'`, registering an `onCleanup` handler that restores the previous visibility and executes `close all force` afterwards.
- Builds the test trajectory with `local_test_trajectory`, then runs `local_test_trajan` and `local_test_trajsimil` in sequence, accumulating results.

### trajan() checks

- **correlation_order**: called with an explicit time axis `[0.0 1.0 2.0]`; asserts exactly 2 line objects (one-spin and two-spin orders) via `test_true`, and, if the count guard holds, verifies with `test_close` (tolerances `1e-14`) that the sorted `XData` of the first line matches the supplied time axis.
- **coherence_order**: asserts exactly 5 line objects, corresponding to coherence orders from -2 to +2 for two spins.
- **total_each_spin**: asserts exactly 2 line objects, one trace per spin.
- **local_each_spin**: asserts exactly 2 line objects, one local trace per spin.
- **level_populations**: asserts exactly 4 line objects, one trace per Zeeman energy level.

Each branch opens a figure with `'Visible','off'`, calls `trajan`, inspects line objects with `findobj(gca,'Type','line')`, and closes the figure.

### trajsimil() checks

Uses the trajectory compared against itself (`traj_ref = traj`) for exact similarity references:

- **RSP** (running scalar product): compares the observed score to `sum(conj(traj).*traj_ref,1)` with `test_close` at tolerances `1e-14`.
- **RDN** (running difference norm): expects unit similarity, compared to `ones(1,size(traj,2))` at tolerances `1e-14`.
- **SG-RDN** (sign-grouped difference norm): expects unit similarity under sign grouping, same tolerances.
- **BSG-RDN** (broad-state-grouped difference norm): expects unit similarity under broad state grouping, same tolerances.

### Test trajectory construction

- Spin system: magnet field 14.1, isotopes `{'1H','13C'}`, scalar Zeeman interactions `{1.0,2.0}`, scalar coupling with `inter.coupling.scalar{1,2}=10.0` and `inter.coupling.scalar{2,2}=0.0`, formalism `'sphten-liouv'`, approximation `'none'`, assembled via `test_spin_system`.
- State `rho_a` is built from `unit_state(spin_system)` plus `0.2*state(spin_system,'Lz','1H')`, `0.1*state(spin_system,'Lz','13C')`, `state(spin_system,'Lx','1H')`, `0.5*state(spin_system,'Ly','13C')`, and `0.25*state(spin_system,{'L+','L-'},{1,2})`.
- The trajectory is `[rho_a 2*rho_a 0.5*rho_a]`, a deterministic mixed-order trajectory from physical states.

## Inputs and outputs

**Syntax**:

```matlab
result = test_dynamic_trajectory_frontends()
```

**Outputs**:

- `result` — regression test result with explanatory messages, accumulated from `test_true` and `test_close` assertions across the `trajan()` and `trajsimil()` branches.

The function takes no inputs.

## References

- `trajan` — trajectory plotting front end exercised in correlation-order, coherence-order, total-each-spin, local-each-spin, and level-population modes.
- `trajsimil` — trajectory similarity scoring front end exercised in RSP, RDN, SG-RDN, and BSG-RDN modes.
- `new_test_result`, `test_true`, `test_close` — regression test harness helpers.
- `test_spin_system`, `unit_state`, `state` — spin system and state construction helpers used to build the test trajectory.
