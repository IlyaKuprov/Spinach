# kernel/overloads/@rcv/times.m

- Signature: `C=times(A,B)`

## Purpose

Scale an RCV sparse matrix by a numeric scalar, whether the matrix is the first or second operand.

## Physical / mathematical content

The stored RCV values are multiplied by the scalar; the operation does not change the row and column index arrays.

## Parameters / inputs

- `A`, `B` - exactly one argument is an RCV sparse matrix; the other is a numeric scalar.

## Outputs

- `C` - RCV sparse matrix with scaled values.

## Implementation structure

After checking the operand types and scalar size, the function multiplies the RCV object's `val` array by the scalar and returns that object. A non-scalar or nonnumeric multiplier is rejected.
