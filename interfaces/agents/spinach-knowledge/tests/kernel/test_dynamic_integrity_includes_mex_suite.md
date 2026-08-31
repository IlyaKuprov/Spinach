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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Announce the test target; implemented by `fprintf('TESTING: Dynamic integrity, include, and MEX helper coverage\n')`.
- Lines 21-24: State the integrity/include/MEX target of the test; implemented by `result=new_test_result('kernel/dynamic_integrity_includes_mex', 'Dynamic integrity, include, and MEX helper coverage', 'low-feasibility include scripts and integrity uti…`.
- Lines 26-27: Locate canonical Spinach subtrees; implemented by `spinach_root=fileparts(fileparts(fileparts(mfilename('fullpath'))))`.
- Lines 32-33: Exercise host overrides and GPU guard includes; implemented by `result=local_test_autoexec(result,fullfile(includes_dir,'autoexec.m'))`.
- Lines 37-38: Start a temporary pool for the pool-dependent include paths; implemented by `pool_created=local_start_pool_if_needed()`.
- Lines 41-43: Exercise parallel profiler and Redfield include paths; implemented by `result=local_test_parallel_profiler(result,fullfile(includes_dir,'parallel_profiler_start.m'), fullfile(includes_dir,'parallel_profiler_report.m'))`.
- Lines 48-49: Exercise integrity utilities using read-only and temporary fixtures; implemented by `result=local_test_existentials(result,fullfile(integrity_dir,'existentials.m'))`.
- Lines 57-58: Exercise MEX compiler helper without building in the real tree; implemented by `result=local_test_compile_mex(result,fullfile(mex_dir,'compile_mex.m'))`.
- Lines 60-61: Release the temporary pool before returning from the test; implemented by `clear('pool_cleanup')`.

### Key state/data transformations

- Lines 22-24: computes `result` using `result=new_test_result('kernel/dynamic_integrity_includes_mex', 'Dynamic integrity, include, and MEX helper coverage', 'low-feasibility include scripts and integrity uti…`.
- Lines 27: computes `spinach_root` using `spinach_root=fileparts(fileparts(fileparts(mfilename('fullpath'))))`.
- Lines 28: computes `includes_dir` using `includes_dir=fullfile(spinach_root,'kernel','includes')`.
- Lines 29: computes `integrity_dir` using `integrity_dir=fullfile(spinach_root,'kernel','integrity')`.
- Lines 30: computes `mex_dir` using `mex_dir=fullfile(spinach_root,'etc','mex')`.
- Lines 38: computes `pool_created` using `pool_created=local_start_pool_if_needed()`.
- Lines 39: computes `pool_cleanup` using `pool_cleanup=onCleanup(@()local_delete_pool(pool_created))`.

### Local helper functions

- Line 66: `local_test_autoexec()` — `function result=local_test_autoexec(result,autoexec_file)`. Preserve environment and graphics defaults touched by autoexec
  - Representative operation: `old_computer=getenv('COMPUTERNAME')`.
  - Representative operation: `old_position=get(groot,'defaultFigurePosition')`.
- Line 123: `local_test_gpu_guard()` — `function result=local_test_gpu_guard(result,start_file,end_file)`. Check GPU removal and restoration around the include pair
  - Representative operation: `spin_system=local_quiet_spin_system()`.
  - Representative operation: `spin_system.sys.enable={'gpu','mex'}`.
- Line 158: `local_test_parallel_profiler()` — `function result=local_test_parallel_profiler(result,start_file,report_file)`. Build a quiet spin system without the detailed dafuq profiler
  - Representative operation: `spin_system=local_quiet_spin_system()`.
  - Representative operation: `spin_system.sys.enable={}`.
- Line 189: `local_test_redfield_serial()` — `function result=local_test_redfield_serial(result,serial_file)`. Build a one-dimensional Redfield fixture with an analytical integral
  - Representative operation: `[spin_system,Q,L0,R,expected]=local_redfield_fixture()`.
  - Representative operation: `run(serial_file)`.
- Line 203: `local_test_redfield_async()` — `function result=local_test_redfield_async(result,async_file)`. Build a one-dimensional Redfield fixture with an analytical integral
  - Representative operation: `[spin_system,Q,L0,R,expected]=local_redfield_fixture()`.
  - Representative operation: `run(async_file)`.
- Line 217: `local_test_direct_include_dispatch()` — `function result=local_test_direct_include_dispatch(result)`. Preserve host and graphics defaults touched by direct include calls
  - Representative operation: `old_computer=getenv('COMPUTERNAME')`.
  - Representative operation: `old_position=get(groot,'defaultFigurePosition')`.
- Line 271: `local_test_existentials()` — `function result=local_test_existentials(result,existentials_file)`. Keep a direct reference to the expected production file
  - Representative operation: `result=test_true(result,'existentials canonical path', strcmp(which('existentials'),existentials_file), 'the read-only existential check must resolve to the canonical Sp…`.
  - Representative operation: `strcmp(which('existentials'),existentials_file), 'the read-only existential check must resolve to the canonical Spinach integrity file')`.
- Line 287: `local_test_exorcise_patrol()` — `function result=local_test_exorcise_patrol(result,exorcise_file,patrol_file)`. Check exorcise input validation without scanning or touching the Wiki
  - Representative operation: `result=test_true(result,'exorcise canonical path',strcmp(which('exorcise'),exorcise_file), 'exorcise must resolve to the production integrity file before validation is t…`.
  - Representative operation: `'exorcise must resolve to the production integrity file before validation is tested')`.

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
