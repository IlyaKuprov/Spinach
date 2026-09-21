# kernel/plotting/axis_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/axis_1d.m`
- Signature: `[ax,ax_label]=axis_1d(spin_system,parameters)`
- Total lines: 220

## Purpose

Generates axis ticks for plotting 1D spectra. Syntax: [ax,ax_label]=axis_1d(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep -either a one-element array giving the sweep width
- in Hz, or a two-element array giving the spectral
- extents in Hz around the offset.
- parameters.zerofill -the number of points in the NMR spectrum after
- zerofilling and Fourier transform
- parameters.offset -offset of the spectrum centre point relative to
- the magnet frequency, Hz
- parameters.axis_units -a character string with the units in which the
- axis ticks should be returned: 'ppm', 'Gauss',
- 'mT', 'Hz', 'kHz', 'MHz','MHz-labframe', 'GHz',
- 'GHz-labframe', 'gtensor', 'points'
- parameters.spins -the spin involved, e.g. {'1H'}

## Outputs

- axis -a row vector of axis tick values
- ax_label -axis label for displaying on the plot
- Note: magnetic field units use the free electron g-tensor for conversion.

## Implementation structure

- Generates axis ticks for plotting 1D spectra. Syntax:
- [ax,ax_label]=axis_1d(spin_system,parameters)
- parameters.sweep - either a one-element array giving the sweep width
- in Hz, or a two-element array giving the spectral
- extents in Hz around the offset.
- parameters.zerofill -the number of points in the NMR spectrum after
- zerofilling and Fourier transform
- parameters.offset -offset of the spectrum centre point relative to
- the magnet frequency, Hz
- parameters.axis_units -a character string with the units in which the
- axis ticks should be returned: 'ppm', 'Gauss',
- 'mT', 'Hz', 'kHz', 'MHz','MHz-labframe', 'GHz',

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `ft_axis()`, `report()`, `num2str()`, `spin()`, `regexp()`, `isfield()`, `ischar()`, `iscell()`, `ismember()`, `strcmp()`.
