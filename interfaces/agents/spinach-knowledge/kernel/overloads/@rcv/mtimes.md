# kernel/overloads/@rcv/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/mtimes.m`
- Signature: `C=mtimes(A,B)`
- Total lines: 90

## Purpose

Multiplication for RCV sparse matrices. Syntax: C=mtimes(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(A,B)`.
- Lines 25-26: RCV sparse by a scalar; implemented by `if isa(A,'rcv')&&isscalar(B)&&isnumeric(B)`.
- Lines 30-31: Scalar by RCV sparse; implemented by `elseif isa(B,'rcv')&&isscalar(A)&&isnumeric(A)`.
- Lines 35-36: RCV sparse by RCV sparse; implemented by `elseif isa(A,'rcv')&&isa(B,'rcv')`.
- Lines 38-39: Check dimension consistency; implemented by `if A.numCols~=B.numRows`.
- Lines 43-44: Result is Matlab sparse; implemented by `C=sparse(A)*sparse(B)`.
- Lines 48-49: Check dimension consistency; implemented by `if A.numCols~=size(B,1)`.
- Lines 53-54: Result is Matlab sparse; implemented by `C=sparse(A)*B`.
- Lines 56-57: Matlab sparse by RCV sparse; implemented by `elseif issparse(A)&&isa(B,'rcv')`.
- Lines 59-60: Check dimension consistency; implemented by `if size(A,2)~=B.numRows`.
- Lines 64-65: Result is Matlab sparse; implemented by `C=A*sparse(B)`.

### Control flow inferred from the code

- Line 26: conditional branch on `isa(A,'rcv')&&isscalar(B)&&isnumeric(B)`.
- Line 39: conditional branch on `A.numCols~=B.numRows`.
- Line 49: conditional branch on `A.numCols~=size(B,1)`.
- Line 60: conditional branch on `size(A,2)~=B.numRows`.

### Key state/data transformations

- Lines 28: computes `A.val` using `A.val=A.val*B; C=A`.
- Lines 33: computes `B.val` using `B.val=B.val*A; C=B`.
- Lines 44: computes `C` using `C=sparse(A)*sparse(B)`.

### Local helper functions

- Line 72: `grumble()` — `function grumble(A,B)`.
  - Representative operation: `if (~isa(A,'rcv'))&&(~isa(B,'rcv'))`.
  - Representative operation: `error('at least one input must be an RCV sparse matrix.')`.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices

## Implementation structure

- Multiplication for RCV sparse matrices. Syntax:
- C=mtimes(A,B)
- A -left operand
- B -right operand
- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- RCV sparse by RCV sparse
- Check dimension consistency
- Result is Matlab sparse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `issparse()`.
