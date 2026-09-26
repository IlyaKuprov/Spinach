# kernel/overloads/@polyadic/suffix.m

- Signature: `p=suffix(p,a)`

## Purpose

Adds suffix matrices to a polyadic. Anything the polyadic multiplies will first be multiplied by the suffix matri- ces. Syntax: p=suffix(p,a)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -polyadic object
- a -suffix matrix

## Outputs

- p -polyadic object
- Note: a suffix can be a polyadic itself.

## Implementation structure

- Adds suffix matrices to a polyadic. Anything the polyadic
- multiplies will first be multiplied by the suffix matri-
- ces. Syntax:
- p=suffix(p,a)
- p - polyadic object
- a - suffix matrix
- Note: a suffix can be a polyadic itself.
- Check consistency
- Absorb the suffix
- Multiply the last core
- Check the dimensions
- Update suffix array
