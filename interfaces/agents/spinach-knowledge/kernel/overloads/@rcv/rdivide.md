# kernel/overloads/@rcv/rdivide.m

- Signature: `A=rdivide(A,k)`

## Purpose

Divides an RCV sparse matrix by a numeric scalar.

## Mathematical content

The operation divides each stored value of `A` by `k`, leaving the sparse structure and dimensions unchanged.

## Numerical / algorithmic content

The function requires `A` to be an RCV object and `k` to be a numeric scalar, then replaces `A.val` with `A.val/k`.

## Parameters / inputs

- A -RCV sparse matrix
- k -numeric scalar

## Outputs

- A -RCV sparse matrix

## Implementation structure

- Check that the first argument is RCV and the divisor is a numeric scalar.
- Divide the stored value array by the scalar and return the modified RCV object.
