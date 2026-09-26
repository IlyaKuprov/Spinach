# kernel/plotting/plot_1d.m

- Signature: `plot_1d(spin_system,spectrum,parameters,varargin)`

## Purpose

1D plotting utility. Syntax: plot_1d(spin_system,spectrum,parameters,varargin)

## Physical / mathematical content

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

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
