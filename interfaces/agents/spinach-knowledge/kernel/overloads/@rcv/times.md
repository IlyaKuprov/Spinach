# kernel/overloads/@rcv/times.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/times.m`
- Signature: `C=times(A,B)`
- Total lines: 57

## Purpose

Multiplies an RCV sparse matrix by a numeric scalar, in either operand order. Syntax: C=times(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(A,B)`.
- Lines 24-25: RCV sparse by a scalar; implemented by `if isa(A,'rcv')`.
- Lines 29-30: Scalar by RCV sparse; implemented by `else`.

### Control flow inferred from the code

- Line 25: conditional branch on `isa(A,'rcv')`.

### Key state/data transformations

- Lines 27: computes `A.val` using `A.val=A.val*B; C=A`.
- Lines 32: computes `B.val` using `B.val=B.val*A; C=B`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(A,B)`.
  - Representative operation: `if ~xor(isa(A,'rcv'),isa(B,'rcv'))`.
  - Representative operation: `error('exactly one of the arguments must be an RCV sparse matrix.')`.

## Parameters / inputs

- A,B -an RCV sparse matrix and a numeric
- scalar, in either order

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Multiplies an RCV sparse matrix by a numeric scalar,
- in either operand order. Syntax:
- C=times(A,B)
- A,B -an RCV sparse matrix and a numeric
- scalar, in either order
- C -RCV sparse matrix
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- Consistency enforcement
- They say that the fish that gets away
- looks bigger than it really is.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xor()`, `isscalar()`.
