# kernel/plotting/plot_uf.m

- Signature: `plot_uf(spin_system,spectrum_uf,parameters)`

## Purpose

Plotting utility for ultrafast constant-time 2D pulse sequences. Syntax: plot_uf(spin_system,spectrum_uf,parameters)

## Physical / mathematical content

## Numerical / algorithmic content

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
