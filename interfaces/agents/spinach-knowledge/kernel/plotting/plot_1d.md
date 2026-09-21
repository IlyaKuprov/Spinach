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
