# kernel/overloads/@ttclass/truncate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@ttclass/truncate.m`
- Signature: `ttout=truncate(tt)`
- Total lines: 72

## Purpose

Performs right-to-left SVD recompression for a tensor train. This should not be called directly, use shrink.m instead. Syntax: ttout=truncate(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

The absolute Frobenius tolerance is converted to relative accuracy using the magnitude of the global coefficient, so negative and complex coefficients retain their phase without changing the error budget. Use `shrink` for public compression; it handles zero coefficients before truncation.

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- tt -a tensor train object with tt.ntrains=1
- and orthogonalised left-to-right

## Outputs

- ttout -a tensor train object, orthogonalised
- right-to-left
- Note: approximation tolerance (in Frobenius norm) is read from
- tt.tolerance property.

## Header notes

Use shrink rather than calling this internal stage directly. The input is a single left-orthogonalised train; the output is right-orthogonalised, with the absolute Frobenius tolerance taken from tt.tolerance.
