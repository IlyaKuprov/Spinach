# kernel/plotting/stack_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/stack_2d.m`
- Signature: `stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)`
- Total lines: 226

## Purpose

Stack plotting utility for 2D NMR spectra. Syntax: stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 38-39: Default measure of stack line opacity; implemented by `if ~exist('alpha_fun','var'), alpha_fun=@(x)sqrt(norm(x,2)); end`.
- Lines 41-42: Check consistency; implemented by `grumble(spectrum,parameters,stack_dim)`.
- Lines 44-45: If a complex spectrum is received, plot both components; implemented by `if nnz(imag(spectrum))>0`.
- Lines 50-51: Inform the user; implemented by `report(spin_system,'plotting ')`.
- Lines 53-54: Accommodate homonuclear 2D sequences; implemented by `if isscalar(parameters.spins)`.
- Lines 64-65: Build axes and apply offsets; implemented by `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2))`.
- Lines 68-69: Convert the units; implemented by `switch parameters.axis_units`.
- Lines 100-101: Get Cartesian arrays of points in 3D; implemented by `[X,Y]=meshgrid(axis_f2,axis_f1); Z=transpose(spectrum)`.
- Lines 103-104: Line direction; implemented by `switch stack_dim`.
- Lines 108-109: Build patch lines; implemented by `patch_x=cell([size(spectrum,2) 1]); patch_y=cell([size(spectrum,2) 1])`.
- Lines 113-114: Close up patch lines; implemented by `patch_x{n}=[X(:,n); NaN]`.
- Lines 118-119: Get the transparency; implemented by `alpha(n)=alpha_fun(Z(:,n))`.
- Lines 125-126: Build patch lines; implemented by `patch_x=cell([size(spectrum,1) 1]); patch_y=cell([size(spectrum,1) 1])`.
- Lines 130-131: Close up patch lines; implemented by `patch_x{n}=[X(n,:)'; NaN]`.
- Lines 135-136: Get the transparency; implemented by `alpha(n)=alpha_fun(Z(n,:))`.
- Lines 142-143: Normalise transparency; implemented by `alpha=alpha/max(alpha)`.
- Lines 146-147: Draw patch lines; implemented by `for n=1:numel(patch_z)`.

### Control flow inferred from the code

- Line 39: conditional branch on `~exist('alpha_fun','var'), alpha_fun=@(x)sqrt(norm(x,2)); end`.
- Line 45: conditional branch on `nnz(imag(spectrum))>0`.
- Line 54: conditional branch on `isscalar(parameters.spins)`.
- Line 57: conditional branch on `isscalar(parameters.offset)`.
- Line 60: conditional branch on `isscalar(parameters.sweep)`.
- Line 69: dispatches on `parameters.axis_units`; cases `'ppm'`, `'Gauss'`, `'Hz'`, `'kHz'`, `'MHz'`, `'points'`.
- Line 104: dispatches on `stack_dim`; cases `1`, `2`.
- Line 111: `for` loop over `n=1:size(Z,2)`.
- Line 128: `for` loop over `n=1:size(Z,1)`.
- Line 147: `for` loop over `n=1:numel(patch_z)`.

### Key state/data transformations

- Lines 36: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 55: computes `parameters.spins` using `parameters.spins= [parameters.spins parameters.spins]`.
- Lines 58: computes `parameters.offset` using `parameters.offset=[parameters.offset parameters.offset]`.
- Lines 61: computes `parameters.sweep` using `parameters.sweep= [parameters.sweep parameters.sweep]`.
- Lines 65: computes `axis_f1` using `axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2))`.
- Lines 66: computes `axis_f2` using `axis_f2=ft_axis(parameters.offset(2),parameters.sweep(2),size(spectrum,1))`.
- Lines 73: computes `axis_f1_label` using `axis_f1_label=['F1: ' parameters.spins{1} ' chemical shift / ppm']`.
- Lines 74: computes `axis_f2_label` using `axis_f2_label=['F2: ' parameters.spins{2} ' chemical shift / ppm']`.
- Lines 101: computes `[X,Y]` using `[X,Y]=meshgrid(axis_f2,axis_f1); Z=transpose(spectrum)`.
- Lines 109: computes `patch_x` using `patch_x=cell([size(spectrum,2) 1]); patch_y=cell([size(spectrum,2) 1])`.
- Lines 110: computes `patch_z` using `patch_z=cell([size(spectrum,2) 1]); alpha=zeros([size(spectrum,2) 1])`.
- Lines 114: computes `patch_x{n}` using `patch_x{n}=[X(:,n); NaN]`.
- Lines 115: computes `patch_y{n}` using `patch_y{n}=[Y(:,n); NaN]`.
- Lines 116: computes `patch_z{n}` using `patch_z{n}=[Z(:,n); NaN]`.
- Lines 119: computes `alpha(n)` using `alpha(n)=alpha_fun(Z(:,n))`.
- Lines 143: computes `alpha` using `alpha=alpha/max(alpha)`.
- Lines 144: computes `alpha(alpha<0.01)` using `alpha(alpha<0.01)=0.01`.

### Local helper functions

- Line 166: `defaults()` — `function parameters=defaults(spin_system,parameters)`. Consistency enforcement
  - Representative operation: `if ~isfield(parameters,'offset')`.
  - Representative operation: `report(spin_system,'parameters.offset field not set, assuming zero offsets.')`.
- Line 178: `grumble()` — `function grumble(spectrum,parameters,stack_dim)`.
  - Representative operation: `if (~isnumeric(spectrum))||(~ismatrix(spectrum))`.
  - Representative operation: `error('spectrum must be a matrix.')`.

## Parameters / inputs

- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep -one or two sweep widths, Hz
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.offset -one or two transmitter offsets, Hz
- parameters.axis_units -axis units ('ppm','Hz','Gauss')
- stack_dim -stacking dimension, 1 or 2
- alpha_fun -optional function handle that takes
- a spectral slice and returns the al-
- pha value that regulates stack line
- opacity

## Outputs

- this function updates the current figure

## Implementation structure

- Stack plotting utility for 2D NMR spectra. Syntax:
- stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)
- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins
- parameters.offset - one or two transmitter offsets, Hz
- parameters.axis_units - axis units ('ppm','Hz','Gauss')
- stack_dim - stacking dimension, 1 or 2
- alpha_fun - optional function handle that takes
- a spectral slice and returns the al-
- pha value that regulates stack line

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `exist()`, `grumble()`, `nnz()`, `report()`, `isscalar()`, `ft_axis()`, `spin()`, `transpose()`, `alpha()`, `alpha_fun()`, `patch()`, `set()`, `camorbit()`, `kxlabel()`, `kylabel()`.
