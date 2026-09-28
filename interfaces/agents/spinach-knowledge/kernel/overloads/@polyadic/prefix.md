# kernel/overloads/@polyadic/prefix.m

- Signature: `p=prefix(a,p)`

## Purpose

Adds prefix matrices to a polyadic. Anything the polyadic multiplies will subsequently be multiplied by the prefix matrices. Syntax: p=prefix(a,p)

## Physical / mathematical content

For nonscalar `a`, the matrix is prepended to the polyadic prefix buffer, so subsequent multiplication applies it from the left; a scalar is absorbed into the first core of each term.

## Numerical / algorithmic content

## Parameters / inputs

- `a`: prefix matrix (which may itself be a polyadic)
- `p`: polyadic object

## Outputs

- `p`: polyadic object
- Note: a prefix can be a polyadic itself.

## Implementation structure

- Checks dimensions for a nonscalar prefix matrix and prepends it to `p.prefix`.
- For scalar `a`, multiplies it into the first core of every term.
- a - prefix matrix
- p - polyadic object
- Note: a prefix can be a polyadic itself.
- Check consistency
- Absorb the prefix
- Multiply the first core
- Check the dimensions
- Update prefix array
