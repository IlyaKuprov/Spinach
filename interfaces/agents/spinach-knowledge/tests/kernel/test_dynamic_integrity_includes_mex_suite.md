# tests/kernel/test_dynamic_integrity_includes_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_integrity_includes_mex_suite.m`
- Signature: `result=test_dynamic_integrity_includes_mex_suite()`
- Total lines: 629

## Purpose

Tests difficult dynamic coverage for includes, integrity, and MEX helpers. Syntax: result=test_dynamic_integrity_includes_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_test_autoexec()`, `local_test_gpu_guard()`, `local_test_parallel_profiler()`, `local_test_redfield_serial()`, `local_test_redfield_async()`, `local_test_direct_include_dispatch()`, `local_test_existentials()`, `local_test_exorcise_patrol()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test uses direct include execution, read-only integrity probes, and
- temporary-directory fixtures for mutating integrity and MEX helpers. It
- avoids touching production tests, production code, shipped build outputs,
- and repository state.

## Implementation structure

- Tests difficult dynamic coverage for includes, integrity, and MEX helpers. Syntax:
- result=test_dynamic_integrity_includes_mex_suite()
- result -regression test result with explanatory messages
- The test uses direct include execution, read-only integrity probes, and
- temporary-directory fixtures for mutating integrity and MEX helpers. It
- avoids touching production tests, production code, shipped build outputs,
- and repository state.
- Announce the test target
- State the integrity/include/MEX target of the test
- Locate canonical Spinach subtrees
- Exercise host overrides and GPU guard includes
- Start a temporary pool for the pool-dependent include paths

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `fileparts()`, `mfilename()`, `fullfile()`, `local_test_autoexec()`, `local_test_gpu_guard()`, `local_start_pool_if_needed()`, `onCleanup()`, `local_delete_pool()`, `local_test_parallel_profiler()`, `local_test_redfield_serial()`, `local_test_redfield_async()`, `local_test_direct_include_dispatch()`, `local_test_existentials()`, `local_test_exorcise_patrol()`, `local_test_rearm_sniff_fixture()`.
