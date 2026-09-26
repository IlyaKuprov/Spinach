# kernel/overloads/@polyadic/plus.m

- Signature: `c=plus(a,b)`

## Purpose

Polyadic addition operation. Does not perform the actual additi- on, but instead stores the operands as a sum of unopened Kronec- ker products. Syntax: c=plus(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

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
