# tests/list_tests.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/list_tests.m`
- Signature: `manifest=list_tests(varargin)`
- Total lines: 37

## Purpose

Lists Spinach regression tests. Syntax: manifest=list_tests(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Add the test library to the path; implemented by `root_dir=fileparts(mfilename('fullpath'))`.
- Lines 21-22: Parse options; implemented by `options=test_options(varargin{:})`.
- Lines 25-26: Apply substring filter; implemented by `if ~isempty(options.pattern)`.
- Lines 32-33: Print the list; implemented by `for n=1:numel(manifest)`.

### Control flow inferred from the code

- Line 26: conditional branch on `~isempty(options.pattern)`.
- Line 33: `for` loop over `n=1:numel(manifest)`.

### Key state/data transformations

- Lines 18: computes `root_dir` using `root_dir=fileparts(mfilename('fullpath'))`.
- Lines 22: computes `options` using `options=test_options(varargin{:})`.
- Lines 23: computes `manifest` using `manifest=test_manifest()`.
- Lines 27-28: computes `keep` using `keep=contains({manifest.id},options.pattern)| contains({manifest.name},options.pattern)`.

## Parameters / inputs

- varargin -optional name-value pair 'pattern', string

## Outputs

- manifest -structure array with test identifiers and names

## Implementation structure

- Lists Spinach regression tests. Syntax:
- manifest=list_tests(varargin)
- varargin -optional name-value pair 'pattern', string
- manifest -structure array with test identifiers and names
- Add the test library to the path
- Parse options
- Apply substring filter
- Print the list

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fileparts()`, `mfilename()`, `addpath()`, `fullfile()`, `test_options()`, `test_manifest()`, `contains()`, `manifest()`.
