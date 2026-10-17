# kernel/plotting/crop_2d.m

- Source: [kernel/plotting/crop_2d.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/crop_2d.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=crop_2d.m)
- Signature: `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)`

## Purpose

Crops a two-dimensional spectrum to two ppm ranges while retaining whole sampled points and updating the frequency-axis metadata. No interpolation is performed.

## Inputs and axes

- `spec` — numeric 2-D spectrum; dimension 1 is the F1 row axis and dimension 2 is the F2 column axis.
- `parameters.sweep`, `parameters.offset`, and `parameters.spins` — one value per axis or a scalar duplicated for both axes. Sweep and offset are in Hz; spins identify the nuclei for ppm conversion.
- `crop_ranges` — two-cell array, `{[f1_min f1_max],[f2_min f2_max]}`; each pair must be finite, real, two-element, ascending and inside its corresponding ppm axis.

Each Hz axis is built with `ft_axis(offset,sweep,size(spec,dimension))`, then converted to ppm using that axis's nuclear spin and the spin-system magnetic field. Thus ppm axes may ascend or descend, including for negative-gyromagnetic-ratio nuclei.

## Sample selection and updated parameters

Bounds select grid indices by strict `>` comparisons, not by rounding to the nearest ppm value. On an ascending axis, the first index above the lower bound is the left index and the first index above the upper bound is the right index. On a descending axis, the corresponding last indices above the upper and lower bounds are used. The returned slice is inclusive from left through right. Exact-boundary samples are not treated as interpolated endpoints, so the retained grid points need not coincide exactly with the requested ppm bounds; the chosen crossing can extend the crop past an upper boundary by one sample.

The cropped matrix is `spec(l_bound_f1:r_bound_f1,l_bound_f2:r_bound_f2)`. The new `parameters.zerofill` is the retained point count on each axis. Digital resolution is the original sweep divided by the original matrix dimension; each new sweep is that resolution times its retained point count. Offsets are recentred from the first retained Hz point, with a parity correction for odd point counts, so the updated metadata reproduces the selected points on the `ft_axis` grid.

## Outputs

- `spec` — cropped 2-D matrix, with retained F1 rows and F2 columns.
- `parameters` — updated `zerofill`, `sweep`, and `offset`; other fields are carried through.
