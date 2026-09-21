# kernel/overloads/@rcv/horzcat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/horzcat.m`
- Signature: `A=horzcat(A,B)`
- Total lines: 60

## Purpose

Horizontal concatenation for RCV sparse matrices. Syntax: A=horzcat(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -left RCV sparse matrix
- B -right RCV sparse matrix

## Outputs

- A -RCV sparse matrix

## Implementation structure

- Horizontal concatenation for RCV sparse matrices. Syntax:
- A=horzcat(A,B)
- A -left RCV sparse matrix
- B -right RCV sparse matrix
- A -RCV sparse matrix
- Check consistency
- Align locations
- Shift column indices
- Concatenate RCV arrays
- Update column count
- Consistency enforcement
- The back half of your forties is a cursed age. It's not

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gpuArray()`.
