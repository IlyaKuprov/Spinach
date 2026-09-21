# kernel/plotting/plot_uf.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_uf.m`
- Signature: `plot_uf(spin_system,spectrum_uf,parameters)`
- Total lines: 184

## Purpose

Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax: plot_uf(spin_system,spectrum_uf,parameters)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spectrum_uf -a real matrix containing the 2D UF NMR
- spectrum
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.dims -sample dimension, m
- parameters.deltat -time step for the acquisition gradient, s
- parameters.npoints -number of points in the acquisition
- gradient
- parameters.Ga -amplitude of the acquisition gradient, T/m
- parameters.offset -two transmitter offsets for the
- conventional and uf dimension, Hz
- parameters.axis_units -axis units ('ppm' or 'Hz')
- parameters.offset_uf_cov -offset between chemical shifts of a MQ
- along the F1 dimension of a conventional
- and an UF spectra ('ppm' or 'Hz').
- Output:
- a figure with correct axis ticks axes in the UF and conventional
- dimension

## Implementation structure

- Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax:
- plot_uf(spin_system,spectrum_uf,parameters)
- spectrum_uf -a real matrix containing the 2D UF NMR
- spectrum
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.dims -sample dimension, m
- parameters.deltat -time step for the acquisition gradient, s
- parameters.npoints -number of points in the acquisition
- gradient
- parameters.Ga -amplitude of the acquisition gradient, T/m
- parameters.offset -two transmitter offsets for the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `flipud()`, `kxlabel()`, `kylabel()`, `set()`, `ismatrix()`, `isfield()`, `isscalar()`, `ischar()`, `ismember()`.
