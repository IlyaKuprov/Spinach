# kernel/plotting/slice_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/slice_2d.m`
- Signature: `slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`
- Total lines: 221

## Purpose

Contour plotting utility with non-linear adaptive contour spacing and 1D slice extraction using mouse. Syntax: slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 61-62: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 64-65: Check consistency; implemented by `grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 67-68: Accommodate homonuclear 2D sequences; implemented by `if isscalar(parameters.spins)`.
- Lines 78-79: Do contour plotting; implemented by `loc=get(0,'defaultfigureposition'); subplot(1,3,1)`.
- Lines 85-86: Switch off performance warning; implemented by `warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')`.
- Lines 88-89: Get the spectrum extents; implemented by `spec_min=min(spectrum(:))`.
- Lines 92-93: Enter the slicing loop; implemented by `while true()`.
- Lines 95-96: Record pointer position; implemented by `subplot(1,3,1); [x,y]=ginput(1)`.
- Lines 98-99: Compute axis grids; implemented by `[F1,F2]=ndgrid(f1,f2)`.
- Lines 101-102: Create the interpolant; implemented by `S_int=griddedInterpolant(F1,F2,transpose(S),'spline')`.
- Lines 104-105: Compute the traces; implemented by `trace_f1=S_int(ones(size(f2))*x,f2)`.
- Lines 108-109: Call the 1D plotting routine for F1; implemented by `parameters_f1=parameters`.
- Lines 117-118: Call the 1D plotting routine for F2; implemented by `parameters_f2=parameters`.

### Control flow inferred from the code

- Line 68: conditional branch on `isscalar(parameters.spins)`.
- Line 71: conditional branch on `isscalar(parameters.offset)`.
- Line 74: conditional branch on `isscalar(parameters.sweep)`.
- Line 93: `while` loop over `true()`.
- Line 119: conditional branch on `numel(parameters.spins)==2`.
- Line 122: conditional branch on `numel(parameters.offset)==2`.

### Key state/data transformations

- Lines 62: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 69: computes `parameters.spins` using `parameters.spins= [parameters.spins parameters.spins]`.
- Lines 72: computes `parameters.offset` using `parameters.offset=[parameters.offset parameters.offset]`.
- Lines 75: computes `parameters.sweep` using `parameters.sweep= [parameters.sweep parameters.sweep]`.
- Lines 79: computes `loc` using `loc=get(0,'defaultfigureposition'); subplot(1,3,1)`.
- Lines 81-82: computes `[f2,f1,S]` using `[f2,f1,S]=plot_2d(spin_system,spectrum,parameters, ncont,delta,k,ncol,m,signs)`.
- Lines 89: computes `spec_min` using `spec_min=min(spectrum(:))`.
- Lines 90: computes `spec_max` using `spec_max=max(spectrum(:))`.
- Lines 96: computes `subplot(1,3,1); [x,y]` using `subplot(1,3,1); [x,y]=ginput(1)`.
- Lines 99: computes `[F1,F2]` using `[F1,F2]=ndgrid(f1,f2)`.
- Lines 102: computes `S_int` using `S_int=griddedInterpolant(F1,F2,transpose(S),'spline')`.
- Lines 105: computes `trace_f1` using `trace_f1=S_int(ones(size(f2))*x,f2)`.
- Lines 106: computes `trace_f2` using `trace_f2=S_int(f1,ones(size(f1))*y)`.
- Lines 109: computes `parameters_f1` using `parameters_f1=parameters`.
- Lines 110: computes `parameters_f1.spins` using `parameters_f1.spins=parameters.spins(1)`.
- Lines 111: computes `parameters_f1.offset` using `parameters_f1.offset=parameters.offset(1)`.
- Lines 112: computes `parameters_f1.sweep` using `parameters_f1.sweep=parameters.sweep(1)`.
- Lines 113: computes `parameters_f1.zerofill` using `parameters_f1.zerofill=parameters.zerofill(1)`.

### Local helper functions

- Line 135: `defaults()` — `function parameters=defaults(spin_system,parameters)`. Consistency enforcement
  - Representative operation: `if ~isfield(parameters,'offset')`.
  - Representative operation: `report(spin_system,'parameters.offset field not set, assuming zero offsets.')`.
- Line 147: `grumble()` — `function grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
  - Representative operation: `if (~isnumeric(spectrum))||(~ismatrix(spectrum))`.
  - Representative operation: `error('spectrum must be a matrix.')`.

## Parameters / inputs

- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep -one or two sweep widths, Hz
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.offset -one or two transmitter offsets, Hz
- parameters.zerofill -one or two point counts in F1 and F2
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

- this function creates a figure
- Note: the following functions are used to compute contour levels:
- cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
- where smin and smax are computed from the spectrum matrix.

## Implementation structure

- Contour plotting utility with non-linear adaptive contour spacing and
- 1D slice extraction using mouse. Syntax:
- slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)
- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins
- parameters.offset - one or two transmitter offsets, Hz
- parameters.zerofill - one or two point counts in F1 and F2
- parameters.axis_units - axis units ('ppm','Hz','Gauss')
- ncont -the number of contours, a reasonable value is 20
- delta -minimum and maximum elevation (as a fraction of the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `grumble()`, `isscalar()`, `get()`, `subplot()`, `set()`, `loc()`, `plot_2d()`, `spectrum()`, `true()`, `ginput()`, `griddedInterpolant()`, `transpose()`, `S_int()`, `plot_1d()`, `isfield()`.
