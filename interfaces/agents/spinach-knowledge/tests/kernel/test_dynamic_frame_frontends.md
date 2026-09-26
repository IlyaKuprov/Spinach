# tests/kernel/test_dynamic_frame_frontends.m

- Signature: `result=test_dynamic_frame_frontends()`

## Purpose

Tests frame-transformation and averaging front-end kernels. Syntax: result=test_dynamic_frame_frontends()

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Outputs

- result -regression test result with explanatory messages
- The test exercises carrier(), frqoffset(), rotframe(), average(), and
- orientation() on compact systems and analytically controlled limiting
- cases.

## Implementation structure

- Tests frame-transformation and averaging front-end kernels. Syntax:
- result=test_dynamic_frame_frontends()
- result -regression test result with explanatory messages
- The test exercises carrier(), frqoffset(), rotframe(), average(), and
- orientation() on compact systems and analytically controlled limiting
- cases.
- Announce the test target
- State the dynamic frame target of the test
- Check carrier Hamiltonian operator-type algebra
- Check single-channel and duplicate-channel frequency offsets
- Check the zeroth-order rotating-frame transformation
- Check average-Hamiltonian and Krylov-Bogolyubov branches
