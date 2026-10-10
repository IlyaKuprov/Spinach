# examples/nqr/nutation_nqr_iodine.m

- Signature: `nutation_nqr_iodine()`
- Source: [`examples/nqr/nutation_nqr_iodine.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nqr/nutation_nqr_iodine.m)

## Purpose and model

This example calculates a powder nutation response for a single 127I nucleus in zero applied magnetic field. The quadrupolar interaction is constructed by `eeqq2nqi(560e6,0.01,5/2,[0 0 0])`; these are the source's arguments, with the first argument supplied as 560e6 and no unit annotation for that argument in this file. It uses the spherical-tensor Liouville formalism without approximation. Relaxation is set to damping at 1e5, with zero equilibrium and a temperature setting of 298 (the source does not state a unit for the latter).

## Nutation schedule and acquisition

The calculation uses the spherical powder grid `rep_2ang_200pts_sph`, a 127I channel, an L+ coil state, and Lx and Ly operators. The frequency sweep is 83.5–84.5 MHz as shown by the code's MHz plotting conversion; the transmitter is set to 84.0 MHz, and the RF-power setting is `2*pi*1e5` (the source does not annotate its unit). Ten powder simulations with `@nqr_pa` use pulse durations `5e-7*n` seconds for n=1…10, i.e. 0.5–5 μs.

For each duration, the code demodulates by multiplying by the transmitter-offset phase, constructs a 512-point frequency axis, and plots the imaginary spectrum in MHz. The panels are labelled in 0.5 μs increments and use vertical limits of −3.0059e−5 to 3.0059e−5. The source comment estimates calculation time as seconds; that is a source note, not a runtime measured here. This is a simulated single-spin powder NQR nutation series; no imported measurement, spatial gradient, chirp, or SPEN/DOSY encoding is specified.
