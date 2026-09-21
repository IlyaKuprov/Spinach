# tests/kernel/test_dynamic_trajectory_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_trajectory_frontends.m`
- Signature: `result=test_dynamic_trajectory_frontends()`
- Total lines: 147

## Purpose

Tests trajectory-analysis dynamic front-end kernels. Syntax: result=test_dynamic_trajectory_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `local_test_trajan()`, `local_test_trajsimil()`, `local_test_trajectory()`, `local_cleanup()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test exercises trajan() plotting branches and trajsimil() scoring
- branches on a compact two-spin spherical-tensor trajectory.

## Implementation structure

- Tests trajectory-analysis dynamic front-end kernels. Syntax:
- result=test_dynamic_trajectory_frontends()
- result -regression test result with explanatory messages
- The test exercises trajan() plotting branches and trajsimil() scoring
- branches on a compact two-spin spherical-tensor trajectory.
- Announce the test target
- State the dynamic trajectory target of the test
- Force invisible figures during plotting checks
- Build the trajectory used by plotting and similarity checks
- Check all trajan() property branches
- Check all trajsimil() scoring families
- Check correlation-order analysis with an explicit time axis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `get()`, `set()`, `onCleanup()`, `local_cleanup()`, `local_test_trajectory()`, `local_test_trajan()`, `local_test_trajsimil()`, `figure()`, `trajan()`, `findobj()`, `test_true()`, `test_close()`, `line_obj()`, `close()`, `trajsimil()`.
