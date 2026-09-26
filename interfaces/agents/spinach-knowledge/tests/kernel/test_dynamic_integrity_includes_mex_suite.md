# tests/kernel/test_dynamic_integrity_includes_mex_suite.m

- Signature: `result=test_dynamic_integrity_includes_mex_suite()`

## Purpose

Tests difficult dynamic coverage for includes, integrity, and MEX helpers. Syntax: result=test_dynamic_integrity_includes_mex_suite()

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

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
