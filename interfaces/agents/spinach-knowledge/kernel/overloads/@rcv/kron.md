# kernel/overloads/@rcv/kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/kron.m`
- Signature: `C=kron(A,B)`
- Total lines: 58

## Purpose

Kronecker product between two RCV sparse matrices. Syntax: C=kron(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(A,B)`.
- Lines 24-25: Compute the output dimensions; implemented by `newNumRows=A.numRows*B.numRows`.
- Lines 28-29: Align locations; implemented by `if A.isGPU||B.isGPU`.
- Lines 34-35: Build the Cartesian product of indices and values; implemented by `[ia,ib]=ndgrid(1:length(A.val),1:length(B.val))`.
- Lines 41-42: Assemble the output RCV object; implemented by `C=rcv(newRow,newCol,newVal,newNumRows,newNumCols)`.

### Control flow inferred from the code

- Line 29: conditional branch on `A.isGPU||B.isGPU`.

### Key state/data transformations

- Lines 25: computes `newNumRows` using `newNumRows=A.numRows*B.numRows`.
- Lines 26: computes `newNumCols` using `newNumCols=A.numCols*B.numCols`.
- Lines 30: computes `A` using `A=gpuArray(A)`.
- Lines 31: computes `B` using `B=gpuArray(B)`.
- Lines 35: computes `[ia,ib]` using `[ia,ib]=ndgrid(1:length(A.val),1:length(B.val))`.
- Lines 36: computes `ia` using `ia=ia(:); ib=ib(:)`.
- Lines 37: computes `newRow` using `newRow=(A.row(ia)-1)*B.numRows+B.row(ib)`.
- Lines 38: computes `newCol` using `newCol=(A.col(ia)-1)*B.numCols+B.col(ib)`.
- Lines 39: computes `newVal` using `newVal=(A.val(ia)).*(B.val(ib))`.
- Lines 42: computes `C` using `C=rcv(newRow,newCol,newVal,newNumRows,newNumCols)`.
- Lines 43: computes `C.isGPU` using `C.isGPU=A.isGPU||B.isGPU`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(A,B)`. Along the Yangzi River, apes moan ceaselessly. My boat has passed ten thousand mounts briskly.
  - Representative operation: `if ~isa(A,'rcv')||~isa(B,'rcv')`.
  - Representative operation: `error('both inputs must be RCV sparse matrices.')`.

## Parameters / inputs

- A -left RCV sparse matrix
- B -right RCV sparse matrix

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Kronecker product between two RCV sparse matrices. Syntax:
- C=kron(A,B)
- A -left RCV sparse matrix
- B -right RCV sparse matrix
- C -RCV sparse matrix
- Check consistency
- Compute the output dimensions
- Align locations
- Build the Cartesian product of indices and values
- Assemble the output RCV object
- Consistency enforcement
- Along the Yangzi River, apes moan ceaselessly.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gpuArray()`, `rcv()`.
