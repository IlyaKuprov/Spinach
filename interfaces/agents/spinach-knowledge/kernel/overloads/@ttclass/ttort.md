# kernel/overloads/@ttclass/ttort.m

- Signature: `[tt,lognrm]=ttort(tt,direct)`

## Purpose

Performs TT-orthogonalisation for a tensor train (or for each tensor train in a buffered sum). Syntax: [tt,lognrm]=ttort(tt,direct)

## Physical / mathematical content

Tensor-train orthogonalisation, applied to one train or to each train in a buffered sum.

## Numerical / algorithmic content

The routine sweeps through the cores using QR decompositions. `direct=+1` orthogonalises left-to-right; `direct=-1` orthogonalises right-to-left. When `lognrm` is requested, it also normalises each buffered train and returns the natural logarithm of its norm.

## Parameters / inputs

- `tt` — tensor-train object, possibly with buffered sums.
- `direct=+1` — left-to-right orthogonality.
- `direct=-1` — right-to-left orthogonality.

## Outputs

- `tt` — tensor-train object with all terms in each buffered sum orthogonalised in the requested direction.
- `lognrm` — when requested, vector of natural logarithms of the norms of the buffered trains; use this option if a tensor norm may exceed `realmax()=1.7977e+308`.

## Header notes

Normally, this subroutine should not be called directly.
