# examples/nqr/pure_nqr_iodine.m

- Signature: `pure_nqr_iodine()`
- Source: [`examples/nqr/pure_nqr_iodine.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nqr/pure_nqr_iodine.m)

## Purpose and model

This is a simulated powder NQR spectrum for one 127I nucleus at zero applied field. The quadrupolar interaction is set by `eeqq2nqi(560e6,0.01,5/2,[0 0 0])`; the first argument is 560e6 in the source, which does not annotate its unit. The basis uses the spherical-tensor Liouville formalism without approximation. Damping is set to 1e5, equilibrium to zero, and temperature to 298 (no unit is stated for the temperature setting).

## Acquisition and processing

The powder calculation uses the `rep_2ang_200pts_sph` grid, the 127I channel, an L+ coil state, an Lx pulse operator, and a π/2 pulse angle. The source sets a sweep parameter of 5e8, 512 points, and `axis_units='MHz'`; the numeric sweep assignment is reported as written because the example does not annotate its unit at that line. The FID comes from `powder(...,@hp_acquire,...,'labframe')`, is exponentially apodised with parameter 6, and is Fourier transformed; the plotted spectrum is the imaginary part of the shifted transform. The source comment estimates calculation time as seconds; that is a source note, not a runtime measured here.

The file defines a single-spin powder simulation. The file does not import measured data or specify a spatial model, gradient/chirp schedule, SPEN, ultrafast DOSY, or multiple-quantum selection.
