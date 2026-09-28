# kernel/utilities/killdiag.m

- Signature: `spec=killdiag(spec,brush_dim)`

## Purpose

Zero a brush-width band around the diagonal of a two-dimensional spectrum.

## Physical / mathematical content

This utility masks matrix values; it does not interpret them beyond their placement in a 2D spectrum.

## Numerical / algorithmic content

For each column, the diagonal row is mapped proportionally from the column index to the row dimension. The rows within the requested brush width are set to zero, with the interval clipped at the matrix boundaries. This mapping also handles rectangular matrices.

## Parameters / inputs

- `spec` - numeric matrix representing a 2D spectrum.
- `brush_dim` - positive real integer no larger than either matrix dimension; specifies the number of rows in the band to mask at each column.

## Outputs

- `spec` - the input matrix with the diagonal band zeroed.

## Implementation structure

After validation, the routine visits each column, calculates the corresponding diagonal row and brush interval, clips the interval to valid row indices, and assigns zero to those entries.
