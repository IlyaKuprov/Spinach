# kernel/plotting/zoom_3d.m

- Signature: `[density,ext]=zoom_3d(density,ext,zoom_ranges)`

## Purpose

Zooms a 3D probability density cube to user-specified fractional limits along each axis.

## Parameters / inputs

- `density`: Real, three-dimensional probability density cube, with dimensions ordered `[X Y Z]`.
- `ext`: Real, six-element vector of grid extents in Angstrom, ordered `[xmin xmax ymin ymax zmin zmax]`. Each minimum must be less than its corresponding maximum.
- `zoom_ranges`: Real, six-element vector of fractional zoom limits, ordered `[xmin xmax ymin ymax zmin zmax]`; for example, `[0.3 0.6 0.1 0.2 0.5 0.8]`. Values must lie between `0` and `1`, and each minimum must be less than its corresponding maximum.

## Outputs

- `density`: Extracted probability density subcube, with dimensions ordered `[X Y Z]`.
- `ext`: Updated grid extents in Angstrom, ordered `[xmin xmax ymin ymax zmin zmax]`.

## Numerical / algorithmic content

The function checks the inputs, then constructs axis coordinates using `linspace` between each pair of supplied extents, with one coordinate per density element along that axis. It converts fractional lower limits to indices using `floor` and fractional upper limits using `ceil`, clamping the resulting indices to the corresponding array bounds. It extracts the indexed subcube and sets `ext` to the original axis coordinates at the selected boundary indices.

## Reference

- [Spinach `zoom_3d.m` page](https://spindynamics.org/wiki/index.php?title=zoom_3d.m)