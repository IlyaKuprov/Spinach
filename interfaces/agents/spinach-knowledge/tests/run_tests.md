# tests/run_tests.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/run_tests.m`
- Signature: `results=run_tests(varargin)`
- Total lines: 94

## Purpose

Runs the Spinach regression test suite. Syntax: results=run_tests(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fileparts()`, `mfilename()`, `addpath()`, `fullfile()`, `genpath()`, `test_options()`, `test_manifest()`, `contains()`, `manifest()`, `feval()`, `results()`, `strcmp()`, `nnz()`, `failed()`.
