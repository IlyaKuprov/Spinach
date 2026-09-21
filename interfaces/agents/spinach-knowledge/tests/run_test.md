# tests/run_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/run_test.m`
- Signature: `result=run_test(test_id)`
- Total lines: 24

## Purpose

Runs one Spinach regression test by identifier substring. Syntax: result=run_test(test_id)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Parameters / inputs

- test_id -test identifier or unique substring from list_tests()

## Outputs

- result -single test result structure

## Implementation structure

- Runs one Spinach regression test by identifier substring. Syntax:
- result=run_test(test_id)
- test_id -test identifier or unique substring from list_tests()
- result -single test result structure
- Run one matching test verbosely

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `run_tests()`.
