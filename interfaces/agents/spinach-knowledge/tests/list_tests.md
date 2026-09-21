# tests/list_tests.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/list_tests.m`
- Signature: `manifest=list_tests(varargin)`
- Total lines: 37

## Purpose

Lists Spinach regression tests. Syntax: manifest=list_tests(varargin)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
