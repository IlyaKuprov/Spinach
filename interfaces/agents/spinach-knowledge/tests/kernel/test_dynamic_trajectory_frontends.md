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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Dynamic trajectory analysis front ends\n')`.
- Lines 19-22: State the dynamic trajectory target of the test; implemented by `result=new_test_result('kernel/dynamic_trajectory_frontends', 'Dynamic trajectory analysis front ends', 'trajectory plotting and similarity helpers must expose determini…`.
- Lines 24-25: Force invisible figures during plotting checks; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 29-30: Build the trajectory used by plotting and similarity checks; implemented by `[spin_system,traj]=local_test_trajectory()`.
- Lines 32-33: Check all trajan() property branches; implemented by `result=local_test_trajan(result,spin_system,traj)`.
- Lines 35-36: Check all trajsimil() scoring families; implemented by `result=local_test_trajsimil(result,spin_system,traj)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_trajectory_frontends', 'Dynamic trajectory analysis front ends', 'trajectory plotting and similarity helpers must expose determini…`.
- Lines 25: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 27: computes `cleaner` using `cleaner=onCleanup(@()local_cleanup(old_visibility))`.
- Lines 30: computes `[spin_system,traj]` using `[spin_system,traj]=local_test_trajectory()`.

### Local helper functions

- Line 41: `local_test_trajan()` — `function result=local_test_trajan(result,spin_system,traj)`. Check correlation-order analysis with an explicit time axis
  - Representative operation: `fig=figure('Visible','off')`.
  - Representative operation: `time_axis=[0.0 1.0 2.0]`.
- Line 89: `local_test_trajsimil()` — `function result=local_test_trajsimil(result,spin_system,traj)`. Use identical trajectories for exact similarity references
  - Representative operation: `traj_ref=traj`.
  - Representative operation: `score_obs=trajsimil(spin_system,traj,traj_ref,'RSP')`.
- Line 114: `local_test_trajectory()` — `function [spin_system,traj]=local_test_trajectory()`. Build a two-spin spherical-tensor Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H','13C'}`.
- Line 139: `local_cleanup()` — `function local_cleanup(old_visibility)`. Restore figure state after success or failure
  - Representative operation: `close all force`.
  - Representative operation: `set(groot,'defaultFigureVisible',old_visibility)`.

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
