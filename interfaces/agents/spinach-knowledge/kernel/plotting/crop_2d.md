# kernel/plotting/crop_2d.m

- Signature: `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)`

## Purpose

Crops 2D spectra to user-specified ranges (in ppm), respecting the digital resolution. Syntax: [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spec -2D matrix containing the spectrum
- crop_ranges -cropping bounds, supplied in the following
- format: {[f1_min f1_max],[f2_min f2_max]}
- The following subfields are required in the parameters structure:
- parameters.sweep -one or two sweep widths, Hz
- parameters.spins -cell array with one ot two character
- strings specifying the working spins.
- parameters.offset -one or two transmitter offsets, Hz

## Outputs

- spec -2D matrix containing the cropped spectrum
- parameters -the updated parameters structure; the new
- offset, sweep, and zerofill reproduce the
- retained axis points on the ft_axis grid

## Implementation structure

- Crops 2D spectra to user-specified ranges (in ppm), respecting the
- digital resolution. Syntax:
- [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)
- spec -2D matrix containing the spectrum
- crop_ranges -cropping bounds, supplied in the following
- format: {[f1_min f1_max],[f2_min f2_max]}
- The following subfields are required in the parameters structure:
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins.
- parameters.offset - one or two transmitter offsets, Hz
- spec -2D matrix containing the cropped spectrum
