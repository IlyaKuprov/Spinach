# kernel/overloads/@rcv/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/full.m`
- Signature: `A=full(A)`
- Total lines: 38

## Purpose

Converts an RCV sparse matrix into a full matrix. Syntax: A=full(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Check consistency; implemented by `grumble(A)`.
- Lines 22-23: Delegate to Matlab; implemented by `A=full(sparse(A))`.

### Key state/data transformations

- Lines 23: computes `A` using `A=full(sparse(A))`.

### Local helper functions

- Line 28: `grumble()` — `function grumble(A)`. Whenever you find yourself on the side of the majority, it is time to pause and reflect.
  - Representative operation: `if ~isa(A,'rcv')`.
  - Representative operation: `error('the input must be an RCV sparse matrix.')`.

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -a full Matlab matrix

## Implementation structure

- Converts an RCV sparse matrix into a full matrix. Syntax:
- A=full(A)
- A -an RCV sparse matrix
- A -a full Matlab matrix
- Check consistency
- Delegate to Matlab
- Consistency enforcement
- Whenever you find yourself on the side of the
- majority, it is time to pause and reflect.
- Mark Twain

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
