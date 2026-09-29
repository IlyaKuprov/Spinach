# kernel/plotting/zoom_3d.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/zoom_3d.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=zoom_3d.m)

- Signature: `[density,ext]=zoom_3d(density,ext,zoom_ranges)`

## Inputs

- `density`: numeric real 3-D cube, with dimensions ordered `[X Y Z]`.
- `ext`: six real grid extents in Angstrom, ordered `[xmin xmax ymin ymax zmin zmax]`; each minimum must be below its maximum.
- `zoom_ranges`: six real fractions ordered `[xmin xmax ymin ymax zmin zmax]`; each value is between 0 and 1 and each lower fraction is below its upper fraction. The source example is `[0.3 0.6 0.1 0.2 0.5 0.8]`.

## Cropping behaviour

For each axis, the source builds coordinates with `linspace(ext_min,ext_max,n)`. It selects from `max(1,floor(n*lower_fraction))` through `min(n,ceil(n*upper_fraction))`, using MATLAB's one-based array indices, then returns that subcube. The returned `ext` is replaced by the coordinates at the selected endpoints. This is index cropping: the routine does not interpolate or resample the density cube, draw a plot, or write a file.
