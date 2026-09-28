# kernel/grids/vcell_solidangle.m

- Signature: `S=vcell_solidangle(P,K,xyz)`

## Purpose

Returns the solid angle of each spherical Voronoi cell specified by `K`. The optional knot points `xyz` guide selection of the cell containing each node rather than its complement.

## Inputs

- `P`: 3 x m array of Voronoi-cell vertex coordinates; columns must be finite real unit vectors, within `1e-6` in squared norm.
- `K`: cell array of nonempty lists of finite, positive integer indices into `P`.
- `xyz` (optional): 3 x numel(K) array of knot points; columns must be finite real unit vectors, within `1e-6` in squared norm.

## Output

- `S`: solid angle for each cell in `K`.

## Implementation

The function validates its inputs, then calls `one_vcell_solidangle` for each cell's vertices. When `xyz` is supplied, it passes the corresponding knot point to that call. The calculation performed by `one_vcell_solidangle` is not included in this source.

Source link: https://spindynamics.org/wiki/index.php?title=vcell_solidangle.m