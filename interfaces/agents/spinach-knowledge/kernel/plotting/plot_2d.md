# kernel/plotting/plot_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_2d.m`
- Signature: `[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...`
- Total lines: 262

## Purpose

Contour plotting utility with non-linear adaptive contour spacing. The function is useful for NMR data where small cross-peaks must be adequa- tely contoured next to large diagonal peaks. Syntax: [axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,... parameters,ncont,delta,... k,ncol,m,signs)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 70-71: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 73-74: Check consistency; implemented by `grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 76-77: Inform the user; implemented by `report(spin_system,'plotting ')`.
- Lines 79-80: If a complex spectrum is received, plot both components; implemented by `if nnz(imag(spectrum))>0`.
- Lines 82-83: Recursively plot the real part; implemented by `subplot(1,2,1)`.
- Lines 87-88: Recursively plot the imaginary part; implemented by `subplot(1,2,2); plot_2d(spin_system,imag(spectrum),parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 91-92: Return the transposed plotting matrix; implemented by `spectrum=transpose(spectrum); return`.
- Lines 96-97: Determine data extents and get contour levels; implemented by `smax=max(spectrum,[],'all'); smin=min(spectrum,[],'all')`.
- Lines 102-103: Accommodate homonuclear 2D sequences; implemented by `if isscalar(parameters.spins)`.
- Lines 113-114: Build axes and apply offsets; implemented by `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2))`.
- Lines 117-118: Convert the units; implemented by `switch parameters.axis_units`.
- Lines 149-150: Plot the spectrum; implemented by `spectrum=transpose(spectrum)`.
- Lines 154-155: Invert the axes; implemented by `set(gca,'XDir','reverse','YDir','reverse')`.
- Lines 157-158: Label the axes; implemented by `kxlabel(axis_f2_label); kylabel(axis_f1_label)`.
- Lines 160-161: Colour the contours; implemented by `if any(positive_contours)&&any(negative_contours)`.
- Lines 178-179: Draw the color bar; implemented by `if ~ismember('colorbar',spin_system.sys.disable)`.

### Control flow inferred from the code

- Line 80: conditional branch on `nnz(imag(spectrum))>0`.
- Line 103: conditional branch on `isscalar(parameters.spins)`.
- Line 106: conditional branch on `isscalar(parameters.offset)`.
- Line 109: conditional branch on `isscalar(parameters.sweep)`.
- Line 118: dispatches on `parameters.axis_units`; cases `'ppm'`, `'Gauss'`, `'Hz'`, `'kHz'`, `'MHz'`, `'points'`.
- Line 161: conditional branch on `any(positive_contours)&&any(negative_contours)`.
- Line 179: conditional branch on `~ismember('colorbar',spin_system.sys.disable)`.

### Key state/data transformations

- Lines 71: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 84: computes `[axis_f1,axis_f2]` using `[axis_f1,axis_f2]=plot_2d(spin_system,real(spectrum),parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 92: computes `spectrum` using `spectrum=transpose(spectrum); return`.
- Lines 97: computes `smax` using `smax=max(spectrum,[],'all'); smin=min(spectrum,[],'all')`.
- Lines 98-100: computes `[contours, positive_contours, negative_contours]` using `[contours, positive_contours, negative_contours]=contspacing(smax,smin,delta,k,signs,ncont)`.
- Lines 99-100: computes `positive_contours, negative_contours]` using `positive_contours, negative_contours]=contspacing(smax,smin,delta,k,signs,ncont)`.
- Lines 100: computes `negative_contours]` using `negative_contours]=contspacing(smax,smin,delta,k,signs,ncont)`.
- Lines 104: computes `parameters.spins` using `parameters.spins= [parameters.spins parameters.spins]`.
- Lines 107: computes `parameters.offset` using `parameters.offset=[parameters.offset parameters.offset]`.
- Lines 110: computes `parameters.sweep` using `parameters.sweep= [parameters.sweep parameters.sweep]`.
- Lines 114: computes `axis_f1` using `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2))`.
- Lines 115: computes `axis_f2` using `axis_f2=ft_axis(parameters.offset(2),parameters.sweep(2),size(spectrum,1))`.
- Lines 122: computes `axis_f1_label` using `axis_f1_label=['F1: ' parameters.spins{1} ' chemical shift / ppm']`.
- Lines 123: computes `axis_f2_label` using `axis_f2_label=['F2: ' parameters.spins{2} ' chemical shift / ppm']`.
- Lines 162: computes `plot_range` using `plot_range=max(abs(positive_contours))+max(abs(negative_contours))`.
- Lines 163: computes `nredcont` using `nredcont=ceil(ncol*max(abs(positive_contours))/plot_range)`.
- Lines 164: computes `nbluecont` using `nbluecont=ceil(ncol*max(abs(negative_contours))/plot_range)`.
- Lines 174: computes `colors` using `colors=0.9*(1-[linspace(1,0,nbluecont)' linspace(1,0,nbluecont)' linspace(0,0,nbluecont)'`.

### Local helper functions

- Line 186: `defaults()` — `function parameters=defaults(spin_system,parameters)`. Consistency enforcement
  - Representative operation: `if ~isfield(parameters,'offset')`.
  - Representative operation: `report(spin_system,'parameters.offset field not set, assuming zero offsets.')`.
- Line 198: `grumble()` — `function grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
  - Representative operation: `if (~isnumeric(spectrum))||(~ismatrix(spectrum))`.
  - Representative operation: `error('spectrum must be a matrix.')`.

## Parameters / inputs

- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep -one or two sweep widths, Hz
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.offset -one or two transmitter offsets, Hz
- parameters.axis_units -axis units ('ppm','Hz','Gauss')
- ncont -the number of contours, a reasonable value is 20
- delta -minimum and maximum elevation (as a fraction of the
- total intensity) of the contours above the baseline.
- A good starting value is [0.02 0.2 0.02 0.2]. The
- first pair of numbers refers to the positive conto-
- urs and the second pair to the negative ones.
- k -a coefficient that controls the curvature of the contour
- spacing function: k=1 corresponds to linear spacing and
- k>1 bends the spacing curve to increase the sampling den-
- sity near the baseline. A reasonable value is 2.
- ncol -number of colours in the colour map; around 256 is fine
- m -the curvature of the colour map: m=1 corresponds to a li-
- near colour ramp into the red for positive contours, and
- into the blue for negative contours. A reasonable value
- for high-contrast plotting is 6.
- signs -can be set to 'positive', 'negative' or 'both' -this
- will cause the corresponding contours to be plotted.

## Outputs

- a figure is drawn and the following parameters returned
- axis_f1, axis_f2 -F1 and F2 axis ticks for external
- plotting utilities
- spectrum -2D spectrum array for external
- plotting utilities
- Note: the following functions are used to compute contour levels:
- cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
- where smin and smax are computed from the spectrum matrix.

## Implementation structure

- Contour plotting utility with non-linear adaptive contour spacing. The
- function is useful for NMR data where small cross-peaks must be adequa-
- tely contoured next to large diagonal peaks. Syntax:
- [axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...
- parameters,ncont,delta,...
- k,ncol,m,signs)
- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins
- parameters.offset - one or two transmitter offsets, Hz
- parameters.axis_units - axis units ('ppm','Hz','Gauss')

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `grumble()`, `report()`, `nnz()`, `subplot()`, `ktitle()`, `transpose()`, `contspacing()`, `isscalar()`, `ft_axis()`, `spin()`, `set()`, `kxlabel()`, `kylabel()`, `any()`, `colormap()`.
