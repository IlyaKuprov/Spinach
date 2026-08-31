# kernel/overloads/@cell/minus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@cell/minus.m`
- Signature: `C=minus(A,B)`
- Total lines: 72

## Purpose

Subtracts cell arrays element-by-element. Syntax: A=minus(A,B)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Check consistency; implemented by `grumble(A,B)`.
- Lines 22-23: Decide the topology; implemented by `if iscell(A)&&iscell(B)`.
- Lines 25-26: Subtract cell-by-cell; implemented by `for n=1:numel(A)`.
- Lines 33-34: Subtract from each cell; implemented by `for n=1:numel(A)`.
- Lines 41-42: Subtract from each cell; implemented by `for n=1:numel(B)`.
- Lines 49-50: Complain and bomb out; implemented by `error('unsupported argument combination.')`.

### Control flow inferred from the code

- Line 23: conditional branch on `iscell(A)&&iscell(B)`.
- Line 26: `for` loop over `n=1:numel(A)`.
- Line 34: `for` loop over `n=1:numel(A)`.
- Line 42: `for` loop over `n=1:numel(B)`.

### Key state/data transformations

- Lines 27: computes `A{n}` using `A{n}=A{n}-B{n}`.
- Lines 29: computes `C` using `C=A`.
- Lines 43: computes `B{n}` using `B{n}=A-B{n}`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(A,B)`.
  - Representative operation: `if (~iscell(A))&&(~isnumeric(A))`.
  - Representative operation: `error('A must be either numeric or a cell array.')`.

## Parameters / inputs

- A,B -cell arrays of identical topology

## Outputs

- A -the resulting cell array

## Implementation structure

- Subtracts cell arrays element-by-element. Syntax:
- A=minus(A,B)
- A,B -cell arrays of identical topology
- A -the resulting cell array
- Check consistency
- Decide the topology
- Subtract cell-by-cell
- Subtract from each cell
- Complain and bomb out
- Consistency enforcement
- I can, therefore I am.
- Simone Weil

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `iscell()`, `isequal()`.
