# kernel/overloads/@polyadic/plus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/plus.m`
- Signature: `c=plus(a,b)`
- Total lines: 81

## Purpose

Polyadic addition operation. Does not perform the actual additi- on, but instead stores the operands as a sum of unopened Kronec- ker products. Syntax: c=plus(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(a,b)`.
- Lines 27-28: Run shortcuts; implemented by `if nnz(a)==0, c=simplify(b); return, end`.
- Lines 31-32: Possible cases; implemented by `if ~isa(a,'polyadic')`.
- Lines 34-35: Matrix + polyadic; implemented by `if isempty(b.prefix)&&isempty(b.suffix)`.
- Lines 43-44: Polyadic + matrix; implemented by `if isempty(a.prefix)&&isempty(a.suffix)`.
- Lines 52-54: Polyadic + polyadic; implemented by `if isempty(a.prefix)&&isempty(a.suffix)&& isempty(b.prefix)&&isempty(b.suffix)`.
- Lines 62-63: Simplify the result; implemented by `c=simplify(c)`.

### Control flow inferred from the code

- Line 28: conditional branch on `nnz(a)==0, c=simplify(b); return, end`.
- Line 29: conditional branch on `nnz(b)==0, c=simplify(a); return, end`.
- Line 32: conditional branch on `~isa(a,'polyadic')`.
- Line 35: conditional branch on `isempty(b.prefix)&&isempty(b.suffix)`.
- Line 44: conditional branch on `isempty(a.prefix)&&isempty(a.suffix)`.
- Line 53: conditional branch on `isempty(a.prefix)&&isempty(a.suffix)&&`.

### Key state/data transformations

- Lines 36: computes `c` using `c=b; c.cores=[c.cores {{a}}]`.

### Local helper functions

- Line 68: `grumble()` — `function grumble(a,b)`. The best revenge is massive success.
  - Representative operation: `[nrows_a,ncols_a]=size(a)`.
  - Representative operation: `[nrows_b,ncols_b]=size(b)`.

## Parameters / inputs

- a,b -polyadic objects

## Outputs

- c -polyadic object
- Note: use this operation sparingly -the additions are simply
- buffered, and all subsequent operations will be slower.

## Implementation structure

- Polyadic addition operation. Does not perform the actual additi-
- on, but instead stores the operands as a sum of unopened Kronec-
- ker products. Syntax:
- c=plus(a,b)
- a,b -polyadic objects
- c -polyadic object
- Note: use this operation sparingly -the additions are simply
- buffered, and all subsequent operations will be slower.
- Check consistency
- Run shortcuts
- Possible cases
- Matrix + polyadic

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nnz()`, `simplify()`, `polyadic()`, `isscalar()`.
