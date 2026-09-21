# kernel/overloads/@rcv/vertcat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/vertcat.m`
- Signature: `A=vertcat(A,B)`
- Total lines: 67

## Purpose

Vertical concatenation for RCV sparse matrices. Syntax: A=vertcat(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -top RCV matrix
- B -bottom RCV matrix

## Outputs

- A -concatenated RCV sparse matrix

## Implementation structure

- Vertical concatenation for RCV sparse matrices. Syntax:
- A=vertcat(A,B)
- A -top RCV matrix
- B -bottom RCV matrix
- A -concatenated RCV sparse matrix
- Check consistency
- Align locations
- Shift row indices
- Concatenate indices
- Update row count in the result
- Consistency enforcement
- Frankly speaking, my dear Karl, I do not like this modern word, which all

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gpuArray()`.
