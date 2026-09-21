# tests/lib/test_options.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_options.m`
- Signature: `options=test_options(varargin)`
- Total lines: 62

## Purpose

Parses name-value options for the Spinach test runner. Syntax: options=test_options(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
