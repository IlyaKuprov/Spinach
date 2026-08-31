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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 20-21: Check consistency; implemented by `grumble(A,B)`.
- Lines 23-24: Align locations; implemented by `if A.isGPU||B.isGPU`.
- Lines 29-30: Shift row indices; implemented by `B.row=B.row+A.numRows`.
- Lines 32-33: Concatenate indices; implemented by `A.row=[A.row; B.row]`.
- Lines 37-38: Update row count in the result; implemented by `A.numRows=A.numRows+B.numRows`.

### Control flow inferred from the code

- Line 24: conditional branch on `A.isGPU||B.isGPU`.

### Key state/data transformations

- Lines 25: computes `A` using `A=gpuArray(A)`.
- Lines 26: computes `B` using `B=gpuArray(B)`.
- Lines 30: computes `B.row` using `B.row=B.row+A.numRows`.
- Lines 33: computes `A.row` using `A.row=[A.row; B.row]`.
- Lines 34: computes `A.col` using `A.col=[A.col; B.col]`.
- Lines 35: computes `A.val` using `A.val=[A.val; B.val]`.
- Lines 38: computes `A.numRows` using `A.numRows=A.numRows+B.numRows`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(A,B)`. Frankly speaking, my dear Karl, I do not like this modern word, which all weaklings use to cloak their feelings when they quarrel with the world
  - Representative operation: `if ~isa(A,'rcv')||~isa(B,'rcv')`.
  - Representative operation: `error('both inputs must be RCV sparse matrices.')`.

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
