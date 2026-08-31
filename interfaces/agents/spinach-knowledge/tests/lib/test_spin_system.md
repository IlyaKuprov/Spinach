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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Apply quiet settings used by regression tests; implemented by `sys.output='hush'`.
- Lines 31-32: Build the Spinach object and basis; implemented by `spin_system=create(sys,inter)`.

### Control flow inferred from the code

- Line 23: conditional branch on `isfield(sys,'disable')`.

### Key state/data transformations

- Lines 22: computes `sys.output` using `sys.output='hush'`.
- Lines 24: computes `sys.disable` using `sys.disable=unique([sys.disable {'hygiene'}])`.
- Lines 28: computes `sys.parallel` using `sys.parallel={'local',1}`.
- Lines 29: computes `sys.parprops` using `sys.parprops={}`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.

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
