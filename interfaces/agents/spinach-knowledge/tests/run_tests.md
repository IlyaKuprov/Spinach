# tests/run_tests.m

- Signature: `results=run_tests(varargin)`

## Purpose

Runs the Spinach regression tests, optionally filtered by a substring of test identifiers or names.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `varargin` — name-value options: `pattern`, `verbose`, and `stop_on_fail`.

## Outputs

- `results` — structure array of test outcomes and messages. If any test fails, the function reports the failures and throws an error rather than returning normally.

## Implementation structure

- Adds the test library and Spinach production directories to the MATLAB path, parses options, and loads the test manifest.
- When `pattern` is non-empty, selects tests whose identifier or name contains it; then runs each selected test and records its status, elapsed time, messages, failures, and error. MATLAB errors raised by a test are recorded as failures.
- Continues through failures unless `stop_on_fail` is enabled; `verbose` prints per-test details. Prints a final summary and throws an error if any tests failed.
