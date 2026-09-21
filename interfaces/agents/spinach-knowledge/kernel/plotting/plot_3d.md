# kernel/plotting/plot_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_3d.m`
- Signature: `plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)`
- Total lines: 238

## Purpose

Volume isosurface plotting utility with non-linear adaptive surface spacing. The function plots the 3D volume and the three projections onto the coordinate planes. Syntax: plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spectrum -a real cube containing the 3D NMR spectrum
- parameters.sweep -three sweep widths, Hz
- parameters.spins -cell array with three character
- strings specifying the working
- spins.
- parameters.offset -three transmitter offsets, Hz
- parameters.axis_units -axis units ('ppm','Hz','Gauss')
- nsurf -the number of surfaces, a reasonable value is 20
- delta -minimum and maximum elevation (as a fraction of the
- total intensity) of the surfaces above the baseline.
- A good starting value is [0.02 0.2 0.02 0.2]. The
- first pair of numbers refers to the positive surfa-
- ces and the second pair to the negative ones.
- k -a coefficient that controls the curvature of the surface
- spacing function: k=1 corresponds to linear spacing and
- k>1 bends the spacing curve to increase the sampling den-
- sity near the baseline. A reasonable value is 2.
- signs -can be set to 'positive', 'negative' or 'both' -this
- will cause the corresponding surfaces to be plotted.

## Outputs

- this function creates a figure
- Note: the following functions are used to compute surface levels:
- cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
- where smin and smax are computed from the spectrum cube.

## Implementation structure

- Volume isosurface plotting utility with non-linear adaptive surface
- spacing. The function plots the 3D volume and the three projections
- onto the coordinate planes. Syntax:
- plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)
- spectrum -a real cube containing the 3D NMR spectrum
- parameters.sweep - three sweep widths, Hz
- parameters.spins - cell array with three character
- strings specifying the working
- spins.
- parameters.offset - three transmitter offsets, Hz
- parameters.axis_units - axis units ('ppm','Hz','Gauss')
- nsurf -the number of surfaces, a reasonable value is 20

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `contspacing()`, `ft_axis()`, `spin()`, `subplot()`, `patch()`, `isosurface()`, `set()`, `box()`, `camlight()`, `lighting()`, `kxlabel()`, `kylabel()`, `kzlabel()`, `plot_2d()`.
