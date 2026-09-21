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
