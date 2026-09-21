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
