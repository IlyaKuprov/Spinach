# kernel/overloads/@polyadic/validate.m

- Signature: `validate(p)`

## Purpose

Checks the internal structure of a polyadic object and throws an error if the object does not meet expectations. Syntax: validate(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Implementation structure

- Checks the internal structure of a polyadic object and throws an
- error if the object does not meet expectations. Syntax:
- validate(p)
- p -a polyadic object
- Check the type
- Check core array, level 1
- Check prefix and suffix arrays, level 1
- Check core array, level 2
- Check prefix and suffix arrays, level 2
- Check core dimensions
- Check prefix and suffix dimensions
- Check the number of terms
