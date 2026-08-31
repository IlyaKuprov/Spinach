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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(A,B)`.
- Lines 24-25: Align locations; implemented by `if A.isGPU||B.isGPU`.
- Lines 30-31: Shift column indices; implemented by `B.col=B.col+A.numCols`.
- Lines 33-34: Concatenate RCV arrays; implemented by `A.row=[A.row; B.row]`.
- Lines 38-39: Update column count; implemented by `A.numCols=A.numCols+B.numCols`.

### Control flow inferred from the code

- Line 25: conditional branch on `A.isGPU||B.isGPU`.

### Key state/data transformations

- Lines 26: computes `A` using `A=gpuArray(A)`.
- Lines 27: computes `B` using `B=gpuArray(B)`.
- Lines 31: computes `B.col` using `B.col=B.col+A.numCols`.
- Lines 34: computes `A.row` using `A.row=[A.row; B.row]`.
- Lines 35: computes `A.col` using `A.col=[A.col; B.col]`.
- Lines 36: computes `A.val` using `A.val=[A.val; B.val]`.
- Lines 39: computes `A.numCols` using `A.numCols=A.numCols+B.numCols`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(A,B)`. The back half of your forties is a cursed age. It's not so much that (as the cliché goes) the policemen look
  - Representative operation: `if ~isa(A,'rcv')||~isa(B,'rcv')`.
  - Representative operation: `error('both inputs must be RCV sparse matrices.')`.

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
