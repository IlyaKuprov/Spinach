# kernel/plotting/crop_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/crop_2d.m`
- Signature: `[spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)`
- Total lines: 137

## Purpose

Crops 2D spectra to user-specified ranges (in ppm), respecting the digital resolution. Syntax: [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(spin_system,spec,parameters,crop_ranges)`.
- Lines 39-40: Accommodate homonuclear 2D sequences; implemented by `if isscalar(parameters.spins)`.
- Lines 50-51: Build axes and apply offsets; implemented by `axis_f1_hz=ft_axis(parameters.offset(1),parameters.sweep(1),size(spec,1))`.
- Lines 54-55: Convert the units; implemented by `axis_f1_ppm=1e6*(2*pi)*axis_f1_hz/(spin(parameters.spins{1})*spin_system.inter.magnet)`.
- Lines 58-59: Find array bounds; implemented by `l_bound_f1=find(axis_f1_ppm>crop_ranges{1}(1),1)`.
- Lines 67-68: Digital resolutions of the input grid; implemented by `dres_f1=parameters.sweep(1)/size(spec,1)`.
- Lines 71-72: Update the point counts; implemented by `parameters.zerofill=[(r_bound_f1-l_bound_f1+1) (r_bound_f2-l_bound_f2+1)]`.
- Lines 74-75: Find the new sweeps; implemented by `parameters.sweep=[parameters.zerofill(1)*dres_f1 parameters.zerofill(2)*dres_f2]`.
- Lines 77-79: Find the new offsets; implemented by `parameters.offset=[axis_f1_hz(l_bound_f1)+parameters.sweep(1)/2-mod(parameters.zerofill(1),2)*dres_f1/2 axis_f2_hz(l_bound_f2)+parameters.sweep(2)/2-mod(parameters.zerof…`.
- Lines 81-82: Cut the spectrum; implemented by `spec=spec(l_bound_f1:r_bound_f1,l_bound_f2:r_bound_f2)`.

### Control flow inferred from the code

- Line 40: conditional branch on `isscalar(parameters.spins)`.
- Line 43: conditional branch on `isscalar(parameters.offset)`.
- Line 46: conditional branch on `isscalar(parameters.sweep)`.
- Line 63: conditional branch on `isempty(l_bound_f1)||isempty(r_bound_f1)||isempty(l_bound_f2)||isempty(r_bound_f2)`.

### Key state/data transformations

- Lines 41: computes `parameters.spins` using `parameters.spins= [parameters.spins parameters.spins]`.
- Lines 44: computes `parameters.offset` using `parameters.offset=[parameters.offset parameters.offset]`.
- Lines 47: computes `parameters.sweep` using `parameters.sweep= [parameters.sweep parameters.sweep]`.
- Lines 51: computes `axis_f1_hz` using `axis_f1_hz=ft_axis(parameters.offset(1),parameters.sweep(1),size(spec,1))`.
- Lines 52: computes `axis_f2_hz` using `axis_f2_hz=ft_axis(parameters.offset(2),parameters.sweep(2),size(spec,2))`.
- Lines 55: computes `axis_f1_ppm` using `axis_f1_ppm=1e6*(2*pi)*axis_f1_hz/(spin(parameters.spins{1})*spin_system.inter.magnet)`.
- Lines 56: computes `axis_f2_ppm` using `axis_f2_ppm=1e6*(2*pi)*axis_f2_hz/(spin(parameters.spins{2})*spin_system.inter.magnet)`.
- Lines 59: computes `l_bound_f1` using `l_bound_f1=find(axis_f1_ppm>crop_ranges{1}(1),1)`.
- Lines 60: computes `r_bound_f1` using `r_bound_f1=find(axis_f1_ppm>crop_ranges{1}(2),1)`.
- Lines 61: computes `l_bound_f2` using `l_bound_f2=find(axis_f2_ppm>crop_ranges{2}(1),1)`.
- Lines 62: computes `r_bound_f2` using `r_bound_f2=find(axis_f2_ppm>crop_ranges{2}(2),1)`.
- Lines 68: computes `dres_f1` using `dres_f1=parameters.sweep(1)/size(spec,1)`.
- Lines 69: computes `dres_f2` using `dres_f2=parameters.sweep(2)/size(spec,2)`.
- Lines 72: computes `parameters.zerofill` using `parameters.zerofill=[(r_bound_f1-l_bound_f1+1) (r_bound_f2-l_bound_f2+1)]`.
- Lines 82: computes `spec` using `spec=spec(l_bound_f1:r_bound_f1,l_bound_f2:r_bound_f2)`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(spin_system,spec,parameters,crop_ranges)`.
  - Representative operation: `if (~isnumeric(spec))||(~ismatrix(spec))`.
  - Representative operation: `error('spec must be a matrix.')`.

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
