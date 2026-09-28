# kernel/overloads/@polyadic/transpose.m

- Signature: `p=transpose(p)`

## Purpose

Returns the ordinary transpose of a polyadic matrix representation; it does not take the complex conjugate.

## Physical / mathematical content

For a product of matrices, transposition reverses factor order: `(AB)^T = B^T A^T`. The polyadic representation applies this to each core and reverses the boundary-factor chains.

## Numerical / algorithmic content

The operation transforms the stored factors directly rather than expanding the represented matrix.

## Implementation structure

- Transposes each matrix in every core.
- Forms the new prefix from the reversed old suffix, transposing each factor.
- Forms the new suffix from the reversed old prefix, transposing each factor.
- Stores the two transformed boundary chains on p.
