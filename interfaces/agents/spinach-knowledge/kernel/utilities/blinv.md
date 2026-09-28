# kernel/utilities/blinv.m

- Signature: `[Lsq,Dsq]=blinv(A)`

## Purpose

Compute Blicharski's relaxation-theory invariants from the interaction matrix `A`, following Equations 20-21 in http://doi.org/10.1515/zna-1972-1012. The source notes that an error and a typo in Equation 21 have been corrected and that the function is insensitive to the trace of `A`.

## Parameters / inputs

- `A` - real 3-by-3 interaction matrix.

## Outputs

- `Lsq` - first-rank invariant, calculated as `(A12-A21)^2+(A13-A31)^2+(A23-A32)^2`.
- `Dsq` - second-rank invariant, calculated as `A11^2+A22^2+A33^2-A11*A22-A11*A33-A22*A33+(3/4)*((A12+A21)^2+(A13+A31)^2+(A23+A32)^2)`.
