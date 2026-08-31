# kernel/overloads/@rcv/gpuArray.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@rcv/gpuArray.m`
- Signature: `obj=gpuArray(obj)`
- Total lines: 49

## Purpose

Transfers an RCV sparse matrix to the GPU. Syntax: obj=gpuArray(obj)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Check consistency; implemented by `grumble(obj)`.
- Lines 22-23: Upload to GPU; implemented by `if ~obj.isGPU`.

### Control flow inferred from the code

- Line 23: conditional branch on `~obj.isGPU`.

### Key state/data transformations

- Lines 24: computes `obj.row` using `obj.row=gpuArray(obj.row)`.
- Lines 25: computes `obj.col` using `obj.col=gpuArray(obj.col)`.
- Lines 26: computes `obj.val` using `obj.val=gpuArray(obj.val)`.
- Lines 27: computes `obj.isGPU` using `obj.isGPU=true`.

### Local helper functions

- Line 33: `grumble()` — `function grumble(obj)`. Then it got worse. The book is very, very good. If someone's going to beat you to the punch with a great
  - Representative operation: `if ~isa(obj,'rcv')`.
  - Representative operation: `error('the input must be an RCV sparse matrix.')`.

## Parameters / inputs

- obj -an RCV sparse matrix

## Outputs

- obj -the same matrix with data stored on GPU

## Implementation structure

- Transfers an RCV sparse matrix to the GPU. Syntax:
- obj=gpuArray(obj)
- obj -an RCV sparse matrix
- obj -the same matrix with data stored on GPU
- Check consistency
- Upload to GPU
- Consistency enforcement
- Then it got worse. The book is very, very good. If
- someone's going to beat you to the punch with a great
- book idea, the least they can do is write something
- crap. Not Andrew. Which shouldn't really come as a
- surprise, since the little bastard is prodigiously

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
