# tests/list_tests.m

- Signature: `manifest=list_tests(varargin)`

## Purpose

Returns and prints the available Spinach regression tests, optionally filtered by a substring.

## Parameters / inputs

- `varargin` - optional name-value pair `pattern`, a substring searched in test identifiers and names.

## Outputs

- `manifest` - structure array of test identifiers and names after filtering.

## Implementation structure

- Adds the test library to the path, parses options with `test_options`, obtains the test entries from `test_manifest`, filters on `id` or `name` when `pattern` is non-empty, and prints each remaining identifier and name as a tab-separated line.
