# kernel/overloads/@polyadic/isreal.m

- Signature: `answer=isreal(p)`

## Purpose

Returns true if the polyadic representation is real. Syntax: answer=isreal(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -true if all numeric data in the object is real

## Implementation structure

- Returns true if the polyadic representation is real. Syntax:
- answer=isreal(p)
- p -a polyadic object
- answer -true if all numeric data in the object is real
- Check the core array
- Check prefix and suffix arrays
- All data is real
- A little inaccuracy sometimes saves a ton of explanation.
- H.H. Munro
