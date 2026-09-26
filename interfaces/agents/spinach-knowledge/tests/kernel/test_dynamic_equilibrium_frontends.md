# tests/kernel/test_dynamic_equilibrium_frontends.m

- Signature: `result=test_dynamic_equilibrium_frontends()`

## Purpose

Tests equilibrium and residual-order dynamic front-end kernels. Syntax: result=test_dynamic_equilibrium_frontends()

## Physical / mathematical content

- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
