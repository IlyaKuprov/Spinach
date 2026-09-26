# kernel/overloads/@polyadic/size.m

- Signature: `varargout=size(p,dim)`

## Purpose

Returns the size of the matrix represented by the polyadic. Syntax: answer=size(p,dim)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object
- dim -dimension whose size is required

## Outputs

- answer -a vector with one or two elements

## Implementation structure

- Returns the size of the matrix represented by the polyadic. Syntax:
- answer=size(p,dim)
- p -a polyadic object
- dim -dimension whose size is required
- answer -a vector with one or two elements
- Check consistency
- Get row dimension
- The leftmost matrix in the prefix
- The cores of the polyadic
- Get column dimension
- The rightmost matrix in the suffix
- Compose the answer
