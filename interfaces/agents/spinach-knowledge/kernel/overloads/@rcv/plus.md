# kernel/overloads/@rcv/plus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/plus.m`
- Signature: `C=plus(A,B)`
- Total lines: 98

## Purpose

Adds things to RCV sparse matrices. Syntax: C=plus(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(A,B)`.
- Lines 25-26: Process the special cases; implemented by `if isa(A,'rcv')&&isscalar(B)&&isnumeric(B)`.
- Lines 28-29: Explain the refusal to add a scalar to an RCV object; implemented by `error('adding a scalar makes a matrix non-sparse; inflate manually first.')`.
- Lines 31-32: Add a scalar to an RCV sparse matrix; implemented by `elseif isa(B,'rcv')&&isscalar(A)&&isnumeric(A)`.
- Lines 37-38: Add two RCV sparse matrices; implemented by `elseif isa(A,'rcv')&&isa(B,'rcv')`.
- Lines 40-41: Check for dimension match; implemented by `if A.numRows~=B.numRows||A.numCols~=B.numCols`.
- Lines 45-46: Align locations; implemented by `if A.isGPU||B.isGPU`.
- Lines 51-52: RCV addition; implemented by `A.row=[A.row; B.row]`.
- Lines 56-57: RCV sparse + Matlab sparse; implemented by `elseif isa(A,'rcv')&&issparse(B)`.
- Lines 59-60: Check for dimension match; implemented by `if A.numRows~=size(B,1)||A.numCols~=size(B,2)`.
- Lines 64-65: Recursive call; implemented by `C=plus(A,rcv(B))`.
- Lines 69-70: Check for dimension match; implemented by `if B.numRows~=size(A,1)||B.numCols~=size(A,2)`.
- Lines 74-75: Recursive call; implemented by `C=plus(rcv(A),B)`.

### Control flow inferred from the code

- Line 26: conditional branch on `isa(A,'rcv')&&isscalar(B)&&isnumeric(B)`.
- Line 41: conditional branch on `A.numRows~=B.numRows||A.numCols~=B.numCols`.
- Line 46: conditional branch on `A.isGPU||B.isGPU`.
- Line 60: conditional branch on `A.numRows~=size(B,1)||A.numCols~=size(B,2)`.
- Line 70: conditional branch on `B.numRows~=size(A,1)||B.numCols~=size(A,2)`.

### Key state/data transformations

- Lines 47: computes `A` using `A=gpuArray(A)`.
- Lines 48: computes `B` using `B=gpuArray(B)`.
- Lines 52: computes `A.row` using `A.row=[A.row; B.row]`.
- Lines 53: computes `A.col` using `A.col=[A.col; B.col]`.
- Lines 54: computes `A.val` using `A.val=[A.val; B.val]; C=A`.
- Lines 65: computes `C` using `C=plus(A,rcv(B))`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(A,B)`.
  - Representative operation: `if (~isa(A,'rcv'))&&(~isa(B,'rcv'))`.
  - Representative operation: `error('at least one input must be an RCV sparse matrix.')`.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -sum A+B, RCV sparse matrix

## Implementation structure

- Adds things to RCV sparse matrices. Syntax:
- C=plus(A,B)
- A -left operand
- B -right operand
- C -sum A+B, RCV sparse matrix
- Check consistency
- Process the special cases
- Explain the refusal to add a scalar to an RCV object
- Add a scalar to an RCV sparse matrix
- Add two RCV sparse matrices
- Check for dimension match
- Align locations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `gpuArray()`, `issparse()`, `rcv()`.
