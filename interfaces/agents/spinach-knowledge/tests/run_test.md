# tests/run_test.m

- Signature: `result=run_test(test_id)`

## Purpose

Runs one regression test selected by a unique substring of its identifier or name.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `test_id` — non-empty row character vector; must match exactly one test identifier or name.

## Outputs

- `result` — result structure for the selected test. If the test fails, `run_tests` throws an error and no result is returned.

## Implementation structure

- Validates `test_id`, adds the test library directory to the MATLAB path, and checks for exactly one substring match in the test manifest.
- Calls `run_tests` with that pattern, verbose output enabled, and stop-on-first-failure enabled.
