# kernel/overloads/@polyadic/polyadic.m

- Signature: `p=polyadic(cores)`

## Purpose

Constructs a matrix as a sum of Kronecker-product terms whose factors remain unopened. For example, `cores={{A,B,C},{D,E}}` represents `kron(A,kron(B,C)) + kron(D,E)`. Multiplicative actions can be computed without explicitly expanding those products; the source comment notes that this can save orders of magnitude in CPU time.

## Physical / mathematical content

The `cores` property stores sums of Kronecker-product terms; `prefix` and `suffix` properties hold matrices applied from the left and right.

## Numerical / algorithmic content

- The constructor checks core consistency, stores `cores`, and validates the resulting object.

## Parameters / inputs

- `cores`: a cell array of cell arrays containing matrix factors whose Kronecker products form the terms of the represented matrix.

## Outputs

- `p`: a polyadic representation of the matrix.
- Note: nested polyadics are permitted -the input matrices may be
- polyadics themselves.

## Implementation structure

- Checks that `cores` is consistent.
- Stores the cell array in `p.cores` and validates the constructed object.
