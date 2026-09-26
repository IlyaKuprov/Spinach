# kernel/overloads/@polyadic/simplify.m

- Signature: `p=simplify(p)`

## Purpose

Simplifies the structure of the polyadic object by reordering buffers, dropping inconsequential terms, and flattening nested polyadics where possible. Syntax: p=simplify(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Parameters / inputs

- p -a polyadic object

## Outputs

- p -a polyadic or a numeric object

## Implementation structure

- Simplifies the structure of the polyadic object by reordering buffers,
- dropping inconsequential terms, and flattening nested polyadics where
- possible. Syntax:
- p=simplify(p)
- p -a polyadic object
- p -a polyadic or a numeric object
- Check consistency
- Get size information
- Flush trivial polyadics into all-zero sparse matrices
- Loop until static
- Default disposition
- Simplify prefixes
