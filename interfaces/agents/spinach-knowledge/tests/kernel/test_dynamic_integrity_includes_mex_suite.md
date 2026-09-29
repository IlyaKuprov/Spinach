# tests/kernel/test_dynamic_integrity_includes_mex_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_integrity_includes_mex_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_integrity_includes_mex_suite.m)

## Purpose

Regression test suite covering "difficult" dynamic coverage for Spinach include scripts, integrity utilities, and the MEX compiler helper. It exercises host-specific parallelisation overrides in `autoexec.m`, GPU guard includes, parallel profiler includes, serial and asynchronous Redfield integral includes, the integrity utilities `existentials`, `exorcise`, `patrol`, `rearm`, `sniff`, and `smack`, and `compile_mex.m`. The suite is designed to avoid touching production tests, production code, shipped build outputs, and repository state: it uses direct include execution, read-only integrity probes, and temporary-directory fixtures for mutating integrity and MEX helpers.

## Behaviour

The main function announces the test target, initialises a test result via `new_test_result` with the name `kernel/dynamic_integrity_includes_mex` and the description "Dynamic integrity, include, and MEX helper coverage", then locates the canonical subtrees `kernel/includes`, `kernel/integrity`, and `etc/mex` relative to the Spinach root (three `fileparts` levels above the test file).

Include-script coverage:

- `local_test_autoexec` runs `kernel/includes/autoexec.m` under controlled `COMPUTERNAME` values, checking that the ALAUNDO host with `sys.enable={'gpu'}` sets `sys.parallel={'processes',32}`, that TALOS with GPU enabled sets `sys.parallel={'processes',12}`, that TALOS without GPU arithmetic leaves `sys.parallel` unset, that a user-supplied `sys.parallel={'processes',4}` is preserved, and that an unknown host leaves `sys.parallel` unset. Environment and graphics root defaults (`defaultFigurePosition`, `defaultFigureWindowStyle`, `defaultFigureMenuBar`, `defaultFigureToolbar`) are saved and restored via `onCleanup`.
- `local_test_gpu_guard` runs `start_disallow_gpu.m` and `end_disallow_gpu.m`, verifying that the start include removes only `gpu` from `spin_system.sys.enable` (keeping unrelated flags such as `mex`) while setting `user_wanted_gpu`, that the end include restores `gpu`, that the pair leaves the enable list unchanged when `gpu` was absent, and that `end_disallow_gpu` errors with a message containing `must be preceded by start_disallow_gpu` when run without the start include.
- `local_test_parallel_profiler` runs `parallel_profiler_start.m`, pauses 0.01 seconds, then runs `parallel_profiler_report.m`, checking that `nbytes` is a finite two-element numeric vector and `walltime` is a finite non-negative scalar. It also source-guards the detailed profiler branch by requiring `parallel.internal.profiling.PoolProfiler` in the start file and `parProfiler.drainLog()` in the report file, since that path uses internal MATLAB API code and is not executed.
- `local_test_redfield_serial` and `local_test_redfield_async` run `redfield_integral_serial.m` and `redfield_integral_async.m` against a one-dimensional Redfield fixture with `rlx_onshell=true` and `rlx_shift=0`, comparing the resulting relaxation matrix `R` to an analytical reference with absolute and relative tolerances of `1e-10`. The serial test also checks that the large input cells `Q` and `L0` are cleared afterwards; the asynchronous test checks that the future variable `F` and `ValueStore` handle `store` are cleared.
- `local_test_direct_include_dispatch` re-executes `autoexec`, the GPU guard pair, the parallel profiler pair, and both Redfield includes by script name (rather than `run` on a full path), verifying the same behaviours through direct script dispatch.

Pool handling: `local_start_pool_if_needed` starts a one-worker `Processes` pool via `parpool('Processes',1)` only if no parallel pool exists; an `onCleanup` handle deletes the pool before the test returns, and only if the suite created it.

Integrity-utility coverage:

- `local_test_existentials` verifies that `which('existentials')` resolves to the canonical `kernel/integrity/existentials.m` and captures its output via `evalc`, requiring the text `Running startup checks`.
- `local_test_exorcise_patrol` verifies canonical paths for `exorcise` and `patrol`, checks that `exorcise('bad-mode')` errors with a message containing `online`, and that `patrol` on a nonexistent subject errors with `file list is empty` (with the random number generator state saved and restored). It also source-guards the high-risk full paths: `exorcise.m` must contain `strcmp(mode,'online')` and `checkcode(file_name)`, and `patrol.m` must contain `while hashes_match` and `eval(mfiles(n).name`, because the production patrol owner path is a continuous example runner and is not run in full.
- `local_test_rearm_sniff_fixture` builds a synthetic Spinach-like tree under a temporary directory (copying `rearm.m` and `sniff.m`, writing a dummy function, and saving a placeholder `smells.mat`), shadows the path only inside the fixture, and verifies that `rearm` prints `sniffer rearmed` and populates non-empty `smells` hash records in the fixture `smells.mat`, that `sniff('none')` accepts the freshly rearmed tree and prints `everything smells fine`, and that `sniff('bad-action')` errors mentioning both `none` and `open`. The scratch tree is removed on cleanup.
- `local_test_real_sniff` runs the production `sniff('none')` read-only from the directory containing `smells.mat`, requiring non-empty output.
- `local_test_smack_static` inspects `smack.m` by source only (it is intentionally destructive command-line cleanup code), requiring the presence of `delete(gcp('nocreate'))`, `parcluster('Processes')`, `fclose('all')`, `clear('all')`, and `gpuDeviceCount`.

MEX helper coverage (`local_test_compile_mex`):

- Builds a scratch tree mirroring `etc/mex`, `etc/mex/private`, `kernel/line_shapes`, `kernel/eigenfields`, and `kernel/indexing`, copies `compile_mex.m` into it, and installs a path-local `mex.m` stub under `etc/mex/private` that records all `mex` calls into `mex_calls.mat`.
- Suppresses the `MATLAB:dispatcher:nameConflict` warning (restoring it on cleanup), runs `compile_mex()`, and verifies exactly five recorded compiler calls with the `-R2018a` flag (interleaved-complex API), compiling `lorentzcon.cpp` and `gausscon.cpp` into `kernel/line_shapes`, `cubic_roots.cpp` into `kernel/eigenfields`, and `spsortrows.cpp` and `spunicols.cpp` into `kernel/indexing`, with `-outdir` values matching the corresponding source directories after path-separator normalisation.
- Source-guards the production `compile_mex.m` for the five listed source paths, the `-outdir` specifications `'-outdir',[P '/kernel/line_shapes']`, `'-outdir',[P '/kernel/eigenfields']`, `'-outdir',[P '/kernel/indexing']`, and `'-R2018a'`.
- Checks `compile_mex.m`, `spsortrows.cpp`, and `spunicols.cpp` for the absence of OpenMP markers (`_OPENMP`, `omp_set`, `omp_get`, `omp.h`, `#pragma omp`, `fopenmp`, `libomp`, `gnu_parallel`, `parallel/algorithm`, `OpenMP`), asserting that the sparse indexing MEX compile and source paths do not depend on OpenMP.

The Redfield fixture (`local_redfield_fixture`) constructs a minimal spin system with `spin_system.tols.rlx_integration=1e-6`, `rlx_zero=1e-14`, `prop_chop=1e-14`, `small_matrix=10`, `dense_matrix=0.5`, formalism `sphten-liouv`, `spin_system.rlx.tau_c={1/3}`, a rank-one tensor `Q` with a single nonzero entry at position (2,2), `L0` and `R` as 1-by-1 sparse matrices, and an analytical reference computed as `-weight*(1-exp(-upper_limit))` with `weight=1/3` and `upper_limit=-1.5*(1/rate)*log(1/spin_system.tols.rlx_integration)` for `rate=-1`.

Helper utilities include `local_quiet_spin_system` (a spin system with `output='hush'`, empty enable/disable lists, and `scratch=tempdir`), `local_error_text` (captures an error message from a function handle), `local_write_dummy_function`, `local_restore_autoexec`, and `local_restore_remove` (restores path and directory, clears shadowed functions, and removes the scratch tree).

## Inputs and outputs

```matlab
result = test_dynamic_integrity_includes_mex_suite()
```

**Inputs:** none.

**Outputs:**

- `result` — regression test result structure with explanatory messages, accumulated through `test_true` and `test_close` assertions.

## References

- [Spinach GitHub repository — test_dynamic_integrity_includes_mex_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_integrity_includes_mex_suite.m)
