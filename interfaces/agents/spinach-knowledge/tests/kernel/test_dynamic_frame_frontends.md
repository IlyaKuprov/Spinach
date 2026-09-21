# tests/kernel/test_dynamic_frame_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_frame_frontends.m`
- Signature: `result=test_dynamic_frame_frontends()`
- Total lines: 197

## Purpose

Tests frame-transformation and averaging front-end kernels. Syntax: result=test_dynamic_frame_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file also defines local helper function(s): `local_test_carrier()`, `local_test_frqoffset()`, `local_test_rotframe()`, `local_test_average()`, `local_test_orientation()`, `local_sphten_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_test_carrier()`, `local_test_frqoffset()`, `local_test_rotframe()`, `local_test_average()`, `local_test_orientation()`, `local_sphten_system()`, `carrier()`, `test_close()`, `basefrqs()`, `operator()`, `frqoffset()`, `test_spin_system()`, `assume()`, `rotframe()`, `average()`.
