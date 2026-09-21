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
