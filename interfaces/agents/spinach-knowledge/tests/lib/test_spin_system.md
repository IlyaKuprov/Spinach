# tests/lib/test_spin_system.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/lib/test_spin_system.m`
- Signature: `spin_system=test_spin_system(sys,inter,bas)`
- Total lines: 35

## Purpose

Builds a small quiet Spinach spin system for tests. Syntax: spin_system=test_spin_system(sys,inter,bas)

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- sys -Spinach system specification
- inter -Spinach interaction specification
- bas -Spinach basis specification

## Outputs

- spin_system -Spinach spin system object

## Implementation structure

- Builds a small quiet Spinach spin system for tests. Syntax:
- spin_system=test_spin_system(sys,inter,bas)
- sys -Spinach system specification
- inter -Spinach interaction specification
- bas -Spinach basis specification
- spin_system -Spinach spin system object
- Apply quiet settings used by regression tests
- Build the Spinach object and basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `create()`, `basis()`.
