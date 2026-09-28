# kernel/overloads/@ttclass/truncate.m

- Signature: `ttout=truncate(tt)`

## Purpose

Performs right-to-left SVD recompression for a tensor train. This should not be called directly, use shrink.m instead. Syntax: ttout=truncate(tt)

## Physical / mathematical content

- Right-to-left SVD recompression of a single tensor train.

## Numerical / algorithmic content

The relative approximation accuracy is `tt.tolerance(1,1)/abs(tt.coeff(1,1))/sqrt(tt.ncores)`. The routine sweeps over cores from right to left, truncates each SVD using `frob_chop`, then normalizes the first core and transfers that norm to `ttout.coeff`.

## Parameters / inputs

- `tt` — a single tensor train (`tt.ntrains=1`), orthogonalised left-to-right.
- Approximation tolerance in Frobenius norm is read from `tt.tolerance`.

## Outputs

- `ttout` — the tensor train, orthogonalised right-to-left.

## Header notes

Use shrink rather than calling this internal stage directly. The input is a single left-orthogonalised train; the output is right-orthogonalised, with the absolute Frobenius tolerance taken from tt.tolerance.
