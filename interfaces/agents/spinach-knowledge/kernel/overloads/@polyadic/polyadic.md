# kernel/overloads/@polyadic/polyadic.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/polyadic.m`
- Signature: `p=polyadic(cores)`
- Total lines: 128

## Purpose

Creates an object of a polyadic class. Syntax: p=polyadic(cores) A polyadic is a matrix formed by a Kronecker product, with those products stored unopened. For example, cores={{A,B,C},{D,E}} corresponds to A(x)B(x)C + D(x)E matrix. Any multiplicative acti- on by this matrix may be computed without opening the Kronecker products. This can save orders of magnitude in CPU time.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `iseye()`, `numel()`, `isnumeric()`, `ismatrix()`, `isfloat()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(cores)`.
- Lines 51-52: Store the cores; implemented by `p.cores=cores`.
- Lines 54-55: Validate the result; implemented by `validate(p)`.

### Key state/data transformations

- Lines 52: computes `p.cores` using `p.cores=cores`.

### Local helper functions

- Line 60: `iseye()` — `function n=iseye(pdc)`. Number of elements
  - Representative operation: `n=false()`.
- Line 67: `numel()` — `function answer=numel(p)`. Multiply the dimensions
  - Representative operation: `answer=prod(size(p))`.
- Line 75: `isnumeric()` — `function answer=isnumeric(p)`. Polyadics are matrices
  - Representative operation: `answer=true()`.
- Line 82: `ismatrix()` — `function answer=ismatrix(p)`. Polyadics are floats
  - Representative operation: `answer=true()`.
- Line 89: `isfloat()` — `function answer=isfloat(p)`. Consistency enforcement
  - Representative operation: `answer=true()`.
- Line 100: `grumble()` — `function grumble(cores)`. Sad news arrived from the University of Nottingham as it announced that
  - Representative operation: `if ~iscell(cores)`.
  - Representative operation: `error('cores must be a cell array.')`.

## Parameters / inputs

- cores -a cell array of cell arrays of matrices
- whose Konecker products make up the mat-
- rix of interest.

## Outputs

- p -a polyadic representation of a matrix
- that behaves in many respects like the
- matrix it represents.
- Note: nested polyadics are permitted -the input matrices may be
- polyadics themselves.

## Implementation structure

- Creates an object of a polyadic class. Syntax:
- p=polyadic(cores)
- A polyadic is a matrix formed by a Kronecker product, with those
- products stored unopened. For example,
- cores={{A,B,C},{D,E}}
- corresponds to A(x)B(x)C + D(x)E matrix. Any multiplicative acti-
- on by this matrix may be computed without opening the Kronecker
- products. This can save orders of magnitude in CPU time.
- cores - a cell array of cell arrays of matrices
- whose Konecker products make up the mat-
- rix of interest.
- p - a polyadic representation of a matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `classdef()`, `grumble()`, `validate()`, `iseye()`, `false()`, `true()`, `ismatrix()`, `isfloat()`, `iscell()`.
