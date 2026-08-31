# kernel/overloads/@opium/opium.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@opium/opium.m`
- Signature: `M=opium(dim,coeff)`
- Total lines: 158

## Purpose

Object Pretending It is a Unit Matrix (OPIUM). Syntax: M=opium(dim,coeff)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `nnz()`, `numel()`, `isnumeric()`, `ismatrix()`, `allfinite()`, `iseye()`, `sparse()`, `conj()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(dim,coeff)`.
- Lines 36-37: Store the parameters; implemented by `M.dim=dim; M.coeff=coeff`.

### Key state/data transformations

- Lines 37: computes `M.dim` using `M.dim=dim; M.coeff=coeff`.

### Local helper functions

- Line 42: `nnz()` — `function n=nnz(op)`. Distinguish zero and scaled unit objects
  - Representative operation: `if op.coeff==0`.
  - Representative operation: `n=0`.
- Line 54: `numel()` — `function n=numel(op)`. Numeric property
  - Representative operation: `n=1`.
- Line 61: `isnumeric()` — `function n=isnumeric(op)`. Numeric property
  - Representative operation: `n=true()`.
- Line 68: `ismatrix()` — `function n=ismatrix(op)`. All finite property
  - Representative operation: `n=true()`.
- Line 75: `allfinite()` — `function n=allfinite(op)`. Is actually unit
  - Representative operation: `if isfinite(op.coeff)`.
  - Representative operation: `n=true()`.
- Line 86: `iseye()` — `function n=iseye(op)`. Conversion to sparse
  - Representative operation: `if op.coeff==1`.
  - Representative operation: `n=true()`.
- Line 97: `sparse()` — `function op=sparse(op)`. Return sparse matrix
  - Representative operation: `op=op.coeff*speye(op.dim)`.
- Line 105: `conj()` — `function op=conj(op)`. Conjugate the coefficient
  - Representative operation: `op.coeff=conj(op.coeff)`.

## Parameters / inputs

- dim -dimension of the unit matrix
- coeff -coefficient in front of the unit matrix

## Outputs

- M -an OPIUM representing the specified matrix

## Implementation structure

- Object Pretending It is a Unit Matrix (OPIUM). Syntax:
- M=opium(dim,coeff)
- dim -dimension of the unit matrix
- coeff -coefficient in front of the unit matrix
- M -an OPIUM representing the specified matrix
- Default properties
- Method description
- Constructor function
- Check consistency
- Store the parameters
- Number of non-zeroes
- Distinguish zero and scaled unit objects

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `classdef()`, `grumble()`, `nnz()`, `true()`, `ismatrix()`, `allfinite()`, `false()`, `iseye()`, `speye()`, `conj()`, `ctranspose()`, `gpuArray()`, `gather()`, `isscalar()`.
