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
