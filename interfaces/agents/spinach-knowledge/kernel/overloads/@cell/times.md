# kernel/overloads/@cell/times.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@cell/times.m`
- Signature: `C=times(A,B)`
- Total lines: 84

## Purpose

Multiplies all entries of a cell array by a user-specified scalar or a matching dimension numeric array. Syntax: C=times(A,B)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -a matrix or a cell array thereof
- B -a matrix or a cell array thereof

## Outputs

- C -the resulting cell array
- Note: both arguments cannot be cell arrays.

## Implementation structure

- Multiplies all entries of a cell array by a user-specified
- scalar or a matching dimension numeric array. Syntax:
- C=times(A,B)
- A -a matrix or a cell array thereof
- B -a matrix or a cell array thereof
- C -the resulting cell array
- Note: both arguments cannot be cell arrays.
- Check consistency
- Decide the topology
- Multiply every cell from the left
- Multiply every cell from the right
- Complain and bomb out

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `iscell()`, `isscalar()`, `isequal()`.
