# tests/kernel/test_dynamic_equilibrium_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_equilibrium_frontends.m`
- Signature: `result=test_dynamic_equilibrium_frontends()`
- Total lines: 143

## Purpose

Tests equilibrium and residual-order dynamic front-end kernels. Syntax: result=test_dynamic_equilibrium_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `local_test_thermalize()`, `local_test_steady()`, `local_test_residual()`, `local_liouville_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test exercises thermalize(), steady(), and residual() on compact
- Liouville-space systems with explicit fixed-point references.

## Implementation structure

- Tests equilibrium and residual-order dynamic front-end kernels. Syntax:
- result=test_dynamic_equilibrium_frontends()
- result -regression test result with explanatory messages
- The test exercises thermalize(), steady(), and residual() on compact
- Liouville-space systems with explicit fixed-point references.
- Announce the test target
- State the dynamic equilibrium target of the test
- Check IME and DiBari thermalisation branches
- Check Newton and squaring steady-state solvers
- Check weak residual-order tensor reduction
- Build a one-spin spherical-tensor Liouville-space system
- Thermalise by the inhomogeneous master equation route

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_test_thermalize()`, `local_test_steady()`, `local_test_residual()`, `local_liouville_system()`, `unit_state()`, `state()`, `thermalize()`, `test_close()`, `operator()`, `propagator()`, `rho_ss()`, `steady()`, `test_spin_system()`, `residual()`, `J_after()`.
