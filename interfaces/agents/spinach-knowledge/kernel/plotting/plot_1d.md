# kernel/plotting/plot_1d.m

- Signature: `plot_1d(spin_system,spectrum,parameters,varargin)`

## Purpose

1D plotting utility. Syntax: plot_1d(spin_system,spectrum,parameters,varargin)

## Physical / mathematical content

## Numerical / algorithmic content

- When `parameters.derivative` is enabled, differentiates the spectrum using `fdvec(spectrum,5,1)`.

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

- Applies defaults and validates the inputs, then obtains the plotting axis and label from `axis_1d`.
- If the spectrum is complex, recursively plots its real and imaginary components and adds a legend. Otherwise, optionally differentiates it, plots it with the supplied `varargin`, and sets tight/padded limits, a box, and grid.
- Reverses the x-axis when `parameters.invert_axis` is enabled.

[Source reference](https://spindynamics.org/wiki/index.php?title=plot_1d.m)
