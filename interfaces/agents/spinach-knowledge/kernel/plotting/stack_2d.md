# kernel/plotting/stack_2d.m

- Signature: `stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)`

## Purpose

Stack plotting utility for 2D NMR spectra. Syntax: stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)

## Physical / mathematical content

## Numerical / algorithmic content

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
