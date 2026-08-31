# tests/lib/test_options.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_options.m`
- Signature: `options=test_options(varargin)`
- Total lines: 62

## Purpose

Parses name-value options for the Spinach test runner. Syntax: options=test_options(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Set defaults; implemented by `options.pattern=''`.
- Lines 22-23: Parse name-value pairs; implemented by `if mod(numel(varargin),2)~=0`.

### Control flow inferred from the code

- Line 23: conditional branch on `mod(numel(varargin),2)~=0`.
- Line 26: `for` loop over `n=1:2:numel(varargin)`.
- Line 27: dispatches on `varargin{n}`; cases `'pattern'`, `'verbose'`, `'stop_on_fail'`.
- Line 31: conditional branch on `islogical(varargin{n+1})&&isscalar(varargin{n+1})`.
- Line 45: conditional branch on `islogical(varargin{n+1})&&isscalar(varargin{n+1})`.

### Key state/data transformations

- Lines 18: computes `options.pattern` using `options.pattern=''`.
- Lines 19: computes `options.verbose` using `options.verbose=false`.
- Lines 20: computes `options.stop_on_fail` using `options.stop_on_fail=false`.

## Parameters / inputs

- varargin -name-value option pairs

## Outputs

- options -options structure

## Implementation structure

- Parses name-value options for the Spinach test runner. Syntax:
- options=test_options(varargin)
- varargin -name-value option pairs
- options -options structure
- Set defaults
- Parse name-value pairs

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `islogical()`, `isscalar()`, `ischar()`, `strcmp()`, `isstring()`.
