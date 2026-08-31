# kernel/plotting/plot_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_1d.m`
- Signature: `plot_1d(spin_system,spectrum,parameters,varargin)`
- Total lines: 157

## Purpose

1D plotting utility. Syntax: plot_1d(spin_system,spectrum,parameters,varargin)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Set common defaults; implemented by `parameters=defaults(spin_system,parameters)`.
- Lines 45-46: Check consistency; implemented by `grumble(spin_system,spectrum,parameters)`.
- Lines 48-49: If a complex spectrum is received, plot both components; implemented by `if ~isreal(spectrum)`.
- Lines 51-52: Recursively plot the real and imaginary components; implemented by `plot_1d(spin_system,real(spectrum),parameters,varargin{:}); hold on`.
- Lines 58-59: Get the axis; implemented by `[ax,ax_label]=axis_1d(spin_system,parameters)`.
- Lines 61-62: Compute the derivative if necessary; implemented by `if isfield(parameters,'derivative')&&parameters.derivative`.
- Lines 66-67: Plot the spectrum; implemented by `plot(ax,spectrum,varargin{:})`.
- Lines 71-72: Label the axis; implemented by `kxlabel(ax_label)`.
- Lines 74-75: Invert the axis if necessary; implemented by `if isfield(parameters,'invert_axis')&&parameters.invert_axis`.

### Control flow inferred from the code

- Line 49: conditional branch on `~isreal(spectrum)`.
- Line 62: conditional branch on `isfield(parameters,'derivative')&&parameters.derivative`.
- Line 75: conditional branch on `isfield(parameters,'invert_axis')&&parameters.invert_axis`.

### Key state/data transformations

- Lines 43: computes `parameters` using `parameters=defaults(spin_system,parameters)`.
- Lines 59: computes `[ax,ax_label]` using `[ax,ax_label]=axis_1d(spin_system,parameters)`.
- Lines 63: computes `spectrum` using `spectrum=fdvec(spectrum,5,1)`.

### Local helper functions

- Line 82: `defaults()` — `function parameters=defaults(spin_system,parameters)`.
  - Representative operation: `if (~isfield(parameters,'offset'))&&isscalar(parameters.sweep)`.
  - Representative operation: `report(spin_system,'parameters.offset field not set, assuming zero offsets.')`.
- Line 101: `grumble()` — `function grumble(spin_system,spectrum,parameters)`.
  - Representative operation: `if (~isnumeric(spectrum))`.
  - Representative operation: `error('spectrum must be a numeric array.')`.

## Parameters / inputs

- spectrum a column vector containing the
- spectrum
- parameters.sweep sweep width, Hz
- parameters.spins spin species, e.g. {'1H'}
- parameters.offset transmitter offset, Hz
- parameters.axis_units axis units ('ppm','Gauss',
- 'mT','Hz','kHz','MHz',
- 'MHz-labframe','GHz','GHz-labframe',
- 'gtensor','points')
- parameters.derivative if set to 1, the spectrum is
- differentiated before plotting
- parameters.invert_axis if set to 1, the frequency axis
- is inverted before plotting
- varargin any number of any other para-
- meters; these will be passed
- to Matlab's plot() function

## Outputs

- this function produces a figure

## Implementation structure

- 1D plotting utility. Syntax:
- plot_1d(spin_system,spectrum,parameters,varargin)
- spectrum a column vector containing the
- spectrum
- parameters.sweep sweep width, Hz
- parameters.spins spin species, e.g. {'1H'}
- parameters.offset transmitter offset, Hz
- parameters.axis_units axis units ('ppm','Gauss',
- 'mT','Hz','kHz','MHz',
- 'MHz-labframe','GHz','GHz-labframe',
- 'gtensor','points')
- parameters.derivative if set to 1, the spectrum is

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `grumble()`, `axis_1d()`, `isfield()`, `fdvec()`, `kxlabel()`, `set()`, `isscalar()`, `report()`, `ischar()`, `iscell()`, `ismember()`.
