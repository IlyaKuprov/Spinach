# tests/kernel/test_dynamic_integrity_includes_mex_suite.m

- Signature: `result=test_dynamic_integrity_includes_mex_suite()`

## Purpose

Regression coverage for dynamic include execution, integrity helpers, and MEX helpers.

## Coverage

- Exercises host-specific `autoexec` behavior, GPU guard scripts, direct include dispatch, existential checks, and parallel-profiler includes.
- Tests serial and asynchronous Redfield-integral includes against a one-dimensional fixture.
- Uses read-only integrity probes and temporary-directory fixtures for mutating integrity and MEX helpers.
- The test avoids modifying production tests, production code, shipped build outputs, or repository state.

## Outputs

- `result` — regression-test result with explanatory messages.
