# examples/relaxation_theory/dd_quad_xcorr_1.m

- MATLAB implementation: [examples/relaxation_theory/dd_quad_xcorr_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/dd_quad_xcorr_1.m)

- Signature: `dd_quad_xcorr_1()`.
- Returns: no MATLAB output arguments; displays a simulated one-dimensional spectrum.

## Purpose

A two-spin liquid-state NMR example of Bloch–Redfield–Wangsness (Redfield) relaxation with a dipolar interaction and a quadrupolar interaction. The source notes that Spinach's relaxation module includes their cross-correlations automatically, including dipole–quadrupole cross-correlation; the dipolar interaction is computed from the supplied Cartesian coordinates. The script then simulates and plots an acquisition rather than returning the relaxation superoperator.

## Model and setup

The fixed system uses `sys.magnet=14.1`, isotopes `1H` and `14N`, a scalar-coupling matrix with off-diagonal entries 50, and a quadrupolar tensor with principal values `[1e4 1e4 -2e4]` and zero Euler angles. The two coordinate rows are `[0 0 0]` and `[0 0 1.02]`; the source does not annotate their length units. Relaxation is set to `{'redfield'}`, equilibrium to `'zero'`, retained terms to `'secular'`, and one correlation-time entry to `1e-9`. It uses the complete `sphten-liouv` basis with no approximation, then calls `create` and `basis`.

For acquisition, the observed spin, initial state, and coil are all `1H` (the latter two use `L+`); decoupling is empty and offset is zero. The sweep is 500 Hz, with 128 points and zero-fill to 512 points; the axis is specified in Hz. `liquid(...,@acquire,...,'nmr')` generates the FID, exponential apodisation is applied with parameter 6, and a zero-filled FFT is shifted. The displayed signal is the real part of the spectrum via `plot_1d`.

## Limits

This is one specified two-spin model and one acquisition; it does not vary the couplings, correlation time, or basis, or separately report the cross-correlation contribution. The script's numeric coupling, quadrupolar, coordinate, and correlation-time settings are literals rather than a reusable input interface. Coordinate units are not stated in the source.
