# kernel/plotting/int_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/int_2d.m`
- Signature: `int_2d(spin_system,spectrum,parameters,ncont,...`
- Total lines: 160

## Purpose

Contour plotting utility with non-linear adaptive contour spacing and 2D integration using mouse or an interval file. Syntax: int_2d(spin_system,spectrum,parameters,ncont,... delta,k,ncol,m,signs,filename)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 66-67: Check consistency; implemented by `grumble(filename)`.
- Lines 69-70: Do contour plotting (a detailed grumbler inside); implemented by `[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 72-73: Switch off performance warning; implemented by `warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')`.
- Lines 75-76: Check if the file exists; implemented by `if ~exist(filename,'file')`.
- Lines 78-79: Get the counter going; implemented by `n=1`.
- Lines 81-82: Report to the user; implemented by `report(spin_system,'define the box by clicking on its opposite corners.')`.
- Lines 85-86: Proceed with interactive integration; implemented by `while true()`.
- Lines 88-89: Get mouse input; implemented by `[ranges{n,1},ranges{n,2}]=ginput(2)`.
- Lines 91-92: Write the file; implemented by `save(filename,'ranges'); drawnow`.
- Lines 94-95: Compute axis grids; implemented by `[F1,F2]=ndgrid(f1,f2)`.
- Lines 97-98: Create the interpolant; implemented by `S_int=griddedInterpolant(F1,F2,transpose(S),'spline')`.
- Lines 100-101: Make function handle; implemented by `fhandle=@(x,y)S_int(x,y)`.
- Lines 103-105: Compute the integral; implemented by `I=integral2(fhandle,ranges{n,1}(1),ranges{n,1}(2), ranges{n,2}(1),ranges{n,2}(2),'RelTol',1e-3,'AbsTol',1e-3)`.
- Lines 107-109: Report to the user; implemented by `report(spin_system,['recorded X range: ' num2str(min(ranges{n,1})) ' to ' num2str(max(ranges{n,1}))])`.
- Lines 114-115: Increment counter; implemented by `n=n+1`.
- Lines 121-122: Load ranges from file; implemented by `report(spin_system,'found the ranges file, integrals:')`.
- Lines 125-126: Perform automatic integration; implemented by `for n=1:size(ranges,1)`.
- Lines 141-142: Report to the user; implemented by `report(spin_system,num2str(I))`.

### Control flow inferred from the code

- Line 76: conditional branch on `~exist(filename,'file')`.
- Line 86: `while` loop over `true()`.
- Line 126: `for` loop over `n=1:size(ranges,1)`.

### Key state/data transformations

- Lines 70: computes `[f2,f1,S]` using `[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`.
- Lines 79: computes `n` using `n=1`.
- Lines 89: computes `[ranges{n,1},ranges{n,2}]` using `[ranges{n,1},ranges{n,2}]=ginput(2)`.
- Lines 95: computes `[F1,F2]` using `[F1,F2]=ndgrid(f1,f2)`.
- Lines 98: computes `S_int` using `S_int=griddedInterpolant(F1,F2,transpose(S),'spline')`.
- Lines 101: computes `fhandle` using `fhandle=@(x,y)S_int(x,y)`.
- Lines 104-105: computes `I` using `I=integral2(fhandle,ranges{n,1}(1),ranges{n,1}(2), ranges{n,2}(1),ranges{n,2}(2),'RelTol',1e-3,'AbsTol',1e-3)`.

### Local helper functions

- Line 151: `grumble()` — `function grumble(filename)`. Beware of the person of one book. Thomas Aquinas
  - Representative operation: `if ~ischar(filename)`.
  - Representative operation: `error('filename argument must be a character string.')`.

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
- filename -the name of the interval file. If this file does not
- exist, the integration intervals are saved into this
- file. If it does exist, the function prints the in-
- tegrals over the ontervals found in this file.

## Outputs

- this function creates a figure and writes the interval file
- Note: the following functions are used to compute contour levels:
- cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
- where smin and smax are computed from the spectrum matrix.

## Implementation structure

- Contour plotting utility with non-linear adaptive contour spacing and
- 2D integration using mouse or an interval file. Syntax:
- int_2d(spin_system,spectrum,parameters,ncont,...
- delta,k,ncol,m,signs,filename)
- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins
- parameters.offset - one or two transmitter offsets, Hz
- parameters.axis_units - axis units ('ppm','Hz','Gauss')
- ncont -the number of contours, a reasonable value is 20
- delta -minimum and maximum elevation (as a fraction of the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `plot_2d()`, `exist()`, `report()`, `true()`, `ginput()`, `save()`, `griddedInterpolant()`, `transpose()`, `S_int()`, `integral2()`, `num2str()`, `load()`, `ischar()`.
