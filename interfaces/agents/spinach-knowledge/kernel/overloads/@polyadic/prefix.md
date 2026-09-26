# kernel/overloads/@polyadic/prefix.m

- Signature: `p=prefix(a,p)`

## Purpose

Adds prefix matrices to a polyadic. Anything the polyadic multiplies will subsequently be multiplied by the prefix matrices. Syntax: p=prefix(a,p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- a -prefix matrix
- p -polyadic object

## Outputs

- p -polyadic object
- Note: a prefix can be a polyadic itself.

## Implementation structure

- Adds prefix matrices to a polyadic. Anything the polyadic
- multiplies will subsequently be multiplied by the prefix
- matrices. Syntax:
- p=prefix(a,p)
- a - prefix matrix
- p - polyadic object
- Note: a prefix can be a polyadic itself.
- Check consistency
- Absorb the prefix
- Multiply the first core
- Check the dimensions
- Update prefix array
