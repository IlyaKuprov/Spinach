# kernel/plotting/plot_uf.m

- Signature: `plot_uf(spin_system,spectrum_uf,parameters)`

## Purpose

Plot an ultrafast constant-time 2D NMR spectrum with axes for the UF and conventional dimensions.

## Physical / mathematical content

The conventional-dimension sweep width is `1/(2*Ta)`, where `Ta=parameters.deltat*parameters.npoints`. The UF axis uses the magnetogyric ratio of the first working spin, `k_max=gamma*parameters.Ga*Ta/(2*pi)`, and `constant_c=(-2*(2*parameters.Te))/parameters.dims`. Its resolution is `abs(1/(constant_c*parameters.dims))` and its sweep width is `abs(k_max/constant_c)`.

The source cites *Progress in Nuclear Magnetic Resonance Spectroscopy* **57** (2010), 241 for the constant. See also the [plot_uf.m documentation](https://spindynamics.org/wiki/index.php?title=plot_uf.m).

## Numerical / algorithmic content

The function constructs the conventional and UF frequency axes from the acquisition settings and transmitter offsets. For `ppm` axes, it converts both axes using the respective working spins and `spin_system.inter.magnet`; for `Hz` axes, it retains the frequency values. It then adds `parameters.offset_uf_cov` to the UF axis, plots `flipud(spectrum_uf)` as contours, and reverses both axis directions. The UF dimension of `spectrum_uf` must equal `round(parameters.dims*k_max)`.

## Parameters / inputs

- `spin_system`: supplies `spin_system.inter.magnet` for ppm conversion.
- `spectrum_uf`: real matrix containing the 2D UF NMR spectrum.
- `parameters.spins`: cell array of one or two character strings specifying the working spins; a single spin is used for both dimensions.
- `parameters.dims`: sample dimension, m.
- `parameters.deltat`: acquisition-gradient time step, s.
- `parameters.npoints`: number of points in the acquisition gradient.
- `parameters.nloops`: number of acquisition loops; sets the conventional-axis point count.
- `parameters.Te`: echo time, used to set `t_max=2*parameters.Te`.
- `parameters.Ga`: acquisition-gradient amplitude, T/m.
- `parameters.offset`: two transmitter offsets for the UF and conventional dimensions, Hz.
- `parameters.axis_units`: `'ppm'` or `'Hz'`.
- `parameters.offset_uf_cov`: offset between chemical shifts of a multiple-quantum signal along the F1 dimension of conventional and UF spectra, in the selected axis units.

## Output

A contour figure with axes for the UF and conventional dimensions.