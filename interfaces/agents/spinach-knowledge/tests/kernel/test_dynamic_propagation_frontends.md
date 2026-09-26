# tests/kernel/test_dynamic_propagation_frontends.m

- Signature: `result=test_dynamic_propagation_frontends()`

## Purpose

Tests dynamic propagation front-end kernels on tiny systems. Syntax: result=test_dynamic_propagation_frontends()

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test exercises propagator(), step(), evolution(), krylov(), and
- reduce() against direct finite-dimensional propagation references.

## Implementation structure

- Tests dynamic propagation front-end kernels on tiny systems. Syntax:
- result=test_dynamic_propagation_frontends()
- result -regression test result with explanatory messages
- The test exercises propagator(), step(), evolution(), krylov(), and
- reduce() against direct finite-dimensional propagation references.
- Announce the test target
- State the dynamic propagation target of the test
- Check scaled Taylor propagator and step() branches
- Check evolution() output modes against explicit propagator products
- Check direct krylov() output modes against step() references
- Check reduce() projector invariants and blanket-disable branch
- Build a one-spin Liouville-space system and force Taylor propagation
