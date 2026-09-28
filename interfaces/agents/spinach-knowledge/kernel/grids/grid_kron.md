# kernel/grids/grid_kron.m

- Signature: `[angles,weights]=grid_kron(angles1,weights1,angles2,weights2)`

## Purpose

Constructs the direct product of two spherical grids, tiling one grid using the rotations of the other. Euler angles use the active ZYZ convention and are given in radians as three columns: `[alpha beta gamma]`.

## Parameters / inputs

- `angles1`, `angles2` — Euler-angle matrices for the first and second grids.
- `weights1`, `weights2` — corresponding column vectors of grid weights.

## Outputs

- `angles` — active ZYZ Euler angles of the product grid, in radians.
- `weights` — product-grid weights.

## Implementation

The function converts both angle grids to quaternions, forms every pair of quaternions, multiplies each pair, and converts the results back to Euler angles. It forms the output weights with `kron(weights1,weights2)`. Input checks require real three-column angle matrices, finite real column-vector weights, and matching angle and weight row counts for each grid.

Source reference: <https://spindynamics.org/wiki/index.php?title=grid_kron.m>