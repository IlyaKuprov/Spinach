# tests/list_tests.m

- Signature: `manifest=list_tests(varargin)`

## Purpose

Lists Spinach regression tests. Syntax: manifest=list_tests(varargin)

## Physical / mathematical content

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
