# kernel/overloads/@rcv/rdivide.m

- Signature: `A=rdivide(A,k)`
- Source: [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/rdivide.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rcv/rdivide.m)

## Purpose

Divides an RCV sparse matrix by a numeric scalar.

## Behaviour

The function checks that `A` is an RCV object and `k` is a numeric scalar, then replaces the stored values with `A.val/k`. The row and column coordinate arrays and the matrix dimensions are unchanged, so the result remains in RCV form. The source does not add real, finite, or nonzero restrictions on `k` beyond its numeric-scalar check.

## Inputs and output

- `A` — RCV sparse matrix.
- `k` — numeric scalar.
- Output `A` — the RCV matrix with each stored value divided by `k`.
