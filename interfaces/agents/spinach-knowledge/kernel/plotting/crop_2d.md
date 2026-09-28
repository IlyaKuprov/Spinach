# kernel/plotting/crop_2d.m

- Signature: `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)`

## Purpose

Crops a two-dimensional spectrum to user-specified frequency-axis ranges in ppm while preserving the digital resolution and updating the axis parameters for the retained points.

## Parameters / inputs

- `spin_system` — Spinach spin-system structure, used to convert frequency axes to ppm.
- `spec` — two-dimensional matrix containing the spectrum.
- `parameters.sweep` — one or two sweep widths in Hz.
- `parameters.spins` — cell array containing one or two working-spin isotope names; a single spin is used for both dimensions.
- `parameters.offset` — one or two transmitter offsets in Hz.
- `crop_ranges` — two-element cell array, `{[f1_min f1_max],[f2_min f2_max]}`; each pair gives ascending ppm bounds within its spectrum axis.

## Numerical / algorithmic content

The routine constructs each axis with `ft_axis`, converts it to ppm using the isotope gyromagnetic ratio and the spin-system magnetic field, and selects the array indices bracketing the requested ranges. Bounds outside the available axes are rejected. The returned `parameters.zerofill`, `parameters.sweep`, and `parameters.offset` are recalculated from the retained points and their original digital resolution.

## Outputs

- `spec` — cropped two-dimensional spectrum.
- `parameters` — updated parameters; the new offset, sweep, and zerofill reproduce the retained axis points on the `ft_axis` grid.
