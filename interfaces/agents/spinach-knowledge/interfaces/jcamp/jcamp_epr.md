# interfaces/jcamp/jcamp_epr.m

Exports common native Spinach EMR outputs with
`text=jcamp_epr(spin_system,parameters,signal,kind,info)`.
Field sweeps use the returned `parameters.b_axis` row in tesla and
`mw_freq` in Hz; ENDOR scans use the `n_frq` row in Hz. Both signals are
native rows, explicitly translated into JCAMP column traces without
conjugation or frequency folding. The time/frequency kinds take columns
or 2D [F2,F1] arrays. Time dimensionality comes from npoints or HYSCORE
nsteps; spectral dimensionality comes from zerofill. Actual signal sizes
supply axis lengths, dwell is 1/sweep, and frequency axes use ft_axis.
Scalar sweep/offset values apply to both dimensions.

Numeric signals and named component structures are supported. The function
calls `jcamp_signal`, which constructs EMR SIMULATION blocks and calls the
final writer. Ownership, output filename, CW/PULSE detection, JCAMP method,
physical parameter description, and additional metadata are explicit in
`info`; no detection mode is guessed. Use the explicit-axis function for
DEER delays, arbitrary trajectories, or irregular sampling. No transform,
normalisation, quadrature reconstruction, or interpolation occurs.

[Source](../../../../jcamp/jcamp_epr.m) · [Usage examples](../../../../jcamp/README.md#export-directly-from-spinach-results)
