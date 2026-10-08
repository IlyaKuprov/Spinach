# examples/imaging/dpfgse_signal_suppress.m

## Purpose

Models DPFGSE water suppression for a seven-proton GABA-in-water solution. The example explicitly supplies gradient and shaped-RF settings, runs the spectroscopy imaging sequence, and plots the processed spectrum.

## Spin model

The source sets a 5.9 T magnet, shifts of 3.00, 3.00, 1.88, 1.88, 2.28, 2.28, and 4.80 ppm, and the listed scalar couplings at 7.36 and 7.58 Hz. It uses the `sphten-liouv` formalism with `IK-2`, proximity level 1, and scalar-coupling connectivity. Path tracing and Krylov propagation are disabled.

## Spatial and sequence settings

The spatial sample is configured with `dims=0.30`, `npts=100`, and derivative rule `{'period',3}`. The initial-state phantom is constant, the initial spin state is `Lz`, the receive-coil phantom is uniform, and detection uses `L+`; diffusion is set to zero.

The configured spectrum has offset 800, sweep 1200, 512 acquired points, zero-fill to 2048, and an axis labelled in Hz. The source sets `g_amp=[1e-3 1.5e-3]` and `g_dur=1e-3`, without an inline unit comment for those gradient values. Its water-selective 180-degree pulse uses ten Gaussian-shaped steps, zero phase, frequency entries of 1220, amplitudes scaled by `2*pi*1700`, and equal step durations summing to `2e-3`; the maximum rank is 2. The source does not annotate units for these RF table values or durations, so the literals are retained without assigning units.

## Output and caveats

`imaging` calls `dpfgse_suppress`; the returned FID is exponentially apodised with parameter 6, Fourier transformed with zero-fill 2048, and the real spectrum is plotted. The source estimates minutes of runtime and says a Tesla V100 is faster; its GPU-enable line is commented out, so this file does not itself enable GPU execution. The timing is a source estimate, not a measurement.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/dpfgse_signal_suppress.m)

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
