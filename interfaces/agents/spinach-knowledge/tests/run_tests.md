# tests/run_tests.m

- Signature: `results=run_tests(varargin)`

## Purpose

Runs the Spinach regression test suite. Syntax: results=run_tests(varargin)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- varargin -name-value options: 'pattern', 'verbose', and
- 'stop_on_fail'

## Outputs

- results -structure array with test outcomes and messages

## Implementation structure

- Runs the Spinach regression test suite. Syntax:
- results=run_tests(varargin)
- varargin -name-value options: 'pattern', 'verbose', and
- 'stop_on_fail'
- results -structure array with test outcomes and messages
- Add the test library to the path
- Add the Spinach production directories to the path
- Parse options
- Get the manifest
- Apply substring filter
- Preallocate result array
- Run the tests
