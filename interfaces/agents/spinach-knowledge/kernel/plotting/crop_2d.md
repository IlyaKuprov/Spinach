# kernel/plotting/crop_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/crop_2d.m`
- Signature: `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)`
- Total lines: 153

## Purpose

Crops 2D spectra to user-specified ranges (in ppm), respecting the digital resolution. Syntax: [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `ft_axis()`, `spin()`, `axis_f1_hz()`, `axis_f2_hz()`, `spec()`, `ismatrix()`, `isfield()`, `iscell()`, `ismember()`, `any()`.
