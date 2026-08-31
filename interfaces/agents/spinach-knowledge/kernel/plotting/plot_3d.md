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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 55-56: Check consistency; implemented by `grumble(spectrum,parameters,nsurf,delta,k,signs)`.
- Lines 58-59: Inform the user; implemented by `report(spin_system,'plotting ')`.
- Lines 61-62: Determine data extents and get surface levels; implemented by `xmax=max(spectrum,[],'all'); xmin=min(spectrum,[],'all')`.
- Lines 65-66: Build axes and apply offsets; implemented by `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,1))`.
- Lines 70-71: Convert the units; implemented by `switch parameters.axis_units`.
- Lines 75-76: Account for magnetogyric ratios; implemented by `axis_f1=1000000*(2*pi)*axis_f1/(spin(parameters.spins{1})*spin_system.inter.magnet)`.
- Lines 83-84: Frequency-swept Gauss units; implemented by `axis_f1=10000*(spin_system.inter.magnet-2*pi*axis_f1/spin('E'))`.
- Lines 91-92: No changes necessary; implemented by `axis_f1=1*axis_f1+0; axis_f2=1*axis_f2+0; axis_f3=1*axis_f3+0`.
- Lines 97-98: Complain and bomb out; implemented by `error('unknown axis units.')`.
- Lines 102-103: Plot the 3D spectrum; implemented by `subplot(2,2,1,'replace'); hold('on')`.
- Lines 124-125: Plot F3-F2 projection; implemented by `subplot(2,2,3); proj_params=parameters`.
- Lines 134-135: Plot F2-F1 projection; implemented by `subplot(2,2,4); proj_params=parameters`.
- Lines 144-145: Plot F3-F1 projection; implemented by `subplot(2,2,2); proj_params=parameters`.
- Lines 154-155: Make the figure bigger; implemented by `loc=get(0,'defaultfigureposition')`.

### Control flow inferred from the code

- Line 71: dispatches on `parameters.axis_units`; cases `'ppm'`, `'Gauss'`, `'Hz'`.
- Line 105: `for` loop over `n=surfaces`.

### Key state/data transformations

- Lines 62: computes `xmax` using `xmax=max(spectrum,[],'all'); xmin=min(spectrum,[],'all')`.
- Lines 63: computes `surfaces` using `surfaces=contspacing(xmax,xmin,delta,k,signs,nsurf)`.
- Lines 66: computes `axis_f1` using `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,1))`.
- Lines 67: computes `axis_f2` using `axis_f2=ft_axis(parameters.offset(2),parameters.sweep(2),size(spectrum,2))`.
- Lines 68: computes `axis_f3` using `axis_f3=ft_axis(parameters.offset(3),parameters.sweep(3),size(spectrum,3))`.
- Lines 79: computes `axis_label` using `axis_label='ch. shift / ppm'`.
- Lines 104: computes `[X,Y,Z]` using `[X,Y,Z]=meshgrid(axis_f2,axis_f1,axis_f3)`.
- Lines 106: computes `p` using `p=patch(isosurface(X,Y,Z,spectrum,n))`.
- Lines 125: computes `subplot(2,2,3); proj_params` using `subplot(2,2,3); proj_params=parameters`.
- Lines 126: computes `proj_params.spins` using `proj_params.spins=proj_params.spins([3 2])`.
- Lines 127: computes `proj_params.offset` using `proj_params.offset=proj_params.offset([3 2])`.
- Lines 128: computes `proj_params.sweep` using `proj_params.sweep=proj_params.sweep([3 2])`.
- Lines 129: computes `proj_params.npoints` using `proj_params.npoints=proj_params.npoints([3 2])`.
- Lines 130: computes `proj_params.zerofill` using `proj_params.zerofill=proj_params.zerofill([3 2])`.
- Lines 135: computes `subplot(2,2,4); proj_params` using `subplot(2,2,4); proj_params=parameters`.
- Lines 145: computes `subplot(2,2,2); proj_params` using `subplot(2,2,2); proj_params=parameters`.
- Lines 155: computes `loc` using `loc=get(0,'defaultfigureposition')`.

### Local helper functions

- Line 161: `grumble()` — `function grumble(spectrum,parameters,ncont,delta,k,signs)`.
  - Representative operation: `if (~isnumeric(spectrum))||(~isreal(spectrum))||(ndims(spectrum)~=3)`.
  - Representative operation: `error('spectrum must be a real cube of numbers.')`.

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
