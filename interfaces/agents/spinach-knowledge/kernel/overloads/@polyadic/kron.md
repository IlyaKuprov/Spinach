# kernel/overloads/@polyadic/kron.m

- Signature: `c=kron(a,b)`

## Purpose

Kronecker product function for polyadics. Syntax: c=kron(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

## Parameters / inputs

- a,b -polyadic or numeric objects

## Outputs

- c -polyadic object
- This operation bundles the inputs into a nested polyadic object.

## Implementation structure

- Kronecker product function for polyadics. Syntax:
- c=kron(a,b)
- a,b -polyadic or numeric objects
- c -polyadic object
- This operation bundles the inputs into a nested polyadic object.
- Check consistency
- Put the new term inside the polyadic structure
- Append B to core lists of A
- Prepend A to core lists of B
- Make a nested polyadic
- Simplify
- Consistency enforcement
