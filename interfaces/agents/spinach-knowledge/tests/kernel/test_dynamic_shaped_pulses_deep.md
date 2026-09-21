# tests/kernel/test_dynamic_shaped_pulses_deep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_shaped_pulses_deep.m`
- Signature: `result=test_dynamic_shaped_pulses_deep()`
- Total lines: 140

## Purpose

Tests dynamic shaped-pulse propagation paths. Syntax: result=test_dynamic_shaped_pulses_deep()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouv_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks constant-generator reductions of shaped_pulse_xy
- product-quadrature/method combinations and shaped_pulse_af Fokker-Planck
- propagation paths.

## Implementation structure

- Tests dynamic shaped-pulse propagation paths. Syntax:
- result=test_dynamic_shaped_pulses_deep()
- result -regression test result with explanatory messages
- The test checks constant-generator reductions of shaped_pulse_xy
- product-quadrature/method combinations and shaped_pulse_af Fokker-Planck
- propagation paths.
- Announce the test target
- State the dynamic shaped-pulse target of the test
- Build a one-proton Liouville-space spin system
- Get compact Liouville-space controls
- Define a constant Cartesian RF generator over two slices
- Exercise the Krylov piecewise-constant path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouv_system()`, `operator()`, `state()`, `step()`, `propagator()`, `shaped_pulse_xy()`, `test_close()`, `strcmp()`, `shaped_pulse_af()`, `test_spin_system()`.
