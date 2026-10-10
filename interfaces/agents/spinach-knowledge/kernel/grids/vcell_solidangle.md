# kernel/grids/vcell_solidangle.m

- Signature: `S=vcell_solidangle(P,K,xyz)`, with `xyz` optional

## Purpose

Returns the solid angle for each spherical polygon described by an indexed vertex list. The optional node directions disambiguate the intended cell from its complementary region on the sphere.

## Inputs and output

- `P` is a finite real `3-by-m` array of unit vertex vectors. The norm check allows an absolute squared-norm deviation of `1e-6`.
- `K` is a cell array. Each `K{j}` is a nonempty list of finite positive integer column indices into `P`, within `1:m`; the list order is used as the polygon boundary order.
- Optional `xyz` is a finite real `3-by-numel(K)` array of unit node vectors, one per cell, with the same squared-norm tolerance.
- `S` contains one solid angle per cell, in steradians. Without `xyz`, `cellfun` preserves `K`'s shape; with `xyz`, the result is a column in the order `j=1,...,numel(K)`.

## Calculation

The wrapper selects `P(:,K{j})` for each cell and calls `one_vcell_solidangle`. Without a node direction, that helper triangulates the ordered polygon as a fan from its first vertex. With `xyz(:,j)`, it closes the polygon and triangulates around that node instead, selecting the region containing the node rather than its complement. Each oriented spherical triangle contributes `2*atan2(det(T), 1 + sum(sum(T.*T(:,[2 3 1]),1),2))`; the contributions are summed, so the vertex ordering determines the signed orientation. Input validation checks the stated numeric, shape, finiteness, unit-length, and index constraints; it does not reorder the polygon vertices.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/vcell_solidangle.m)
- [Spinach Wiki: vcell_solidangle.m](https://spindynamics.org/wiki/index.php?title=vcell_solidangle.m)
- The triangulation formula is attributed in the helper source to [DOI: 10.1109/TBME.1983.325207](https://doi.org/10.1109/TBME.1983.325207); see also [one_vcell_solidangle.m on the Spinach Wiki](https://spindynamics.org/wiki/index.php?title=one_vcell_solidangle.m).
