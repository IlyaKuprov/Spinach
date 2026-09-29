# tests/kernel/test_dynamic_iserstep_highorder.m

[View source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_iserstep_highorder.m)

## Purpose

Regression test for the nonlinear high-order `iserstep` branches. It checks zero-step handling, nonlinear generator execution, and agreement of the high-order Lie and RKMK branches against a refined Dormand–Prince 8th-order (DP8) reference on a compact Hilbert-space problem.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression test result via `new_test_result` for `kernel/dynamic_iserstep_highorder`, describing the requirement that high-order Lie and RKMK solvers execute nonlinear generators while preserving density-matrix invariants.
- Builds a one-proton Hilbert-space spin system using `local_hilb_system`, with Pauli operators `Sx`, `Sy`, `Sz` from `pauli(2)`.
- Uses the Hermitian density matrix `rho = [0.7, 0.2+0.1i; 0.2-0.1i, 0.3]`, which has non-zero coherences.
- Defines a mildly nonlinear, non-commuting Hamiltonian field:
  - `H0 = 0.31*Sx + 0.17*Sz`
  - `H1 = -0.23*Sy + 0.05*Sx`
  - `H2 = 0.11*Sx - 0.07*Sz`
  - `Lfun = @(t,state) H0 + sin(4*t)*H1 + 0.2*real(trace(Sz*state))*H2`
  - Time step `dt = 5e-3`.
- Checks the explicit zero-time shortcut in `LG4A`: `iserstep` with `{Lfun, 0, 'LG4A'}` and step `0` must return the input state, compared with absolute and relative tolerances `1e-15`.
- Builds a refined DP8 reference by two consecutive `RKMK-DP8` half steps (`dt/2` each, the second starting at time `dt/2`).
- Exercises the nonlinear high-order branches `LG4A`, `RKMK4`, `RKMK-DP5`, `RKMK-DP8`, and `RKMK-RKF45`, each taking one nonlinear Lie step of size `dt` from time 0 and comparing against the refined reference with per-method tolerances `[5e-8 5e-8 5e-9 5e-10 5e-9]`.
- For each method, checks density-matrix invariants preserved by unitary propagation: trace preservation (`trace(rho_obs)` vs `trace(rho)`) and Hermiticity (`rho_obs` vs `rho_obs'`), both with tolerances `1e-12`.
- All comparisons are accumulated into the test result via `test_close` with explanatory messages.

## Inputs and outputs

```matlab
result = test_dynamic_iserstep_highorder()
```

- **Output**: `result` — regression test result with explanatory messages.
- Takes no inputs.

## References

- Source: [tests/kernel/test_dynamic_iserstep_highorder.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_iserstep_highorder.m) on GitHub.
