# examples/nmr_spen/psyche_strychnine.m

- MATLAB implementation: [examples/nmr_spen/psyche_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/psyche_strychnine.m)

- Signature: `psyche_strychnine()`
- Source: [`examples/nmr_spen/psyche_strychnine.m`](../../../../../examples/nmr_spen/psyche_strychnine.m)

## Model and spatial encoding

The example obtains the proton strychnine model from `strychnine({'1H'})`, shifts the chemical-shift array by −5.0, and sets the magnetic field to 11.7 T. It uses the spherical-tensor Liouville formalism, IK-2 basis approximation, scalar-coupling connectivity, proximal level 1, and an interaction-cutoff tolerance of 2.0. The PSYCHE sequence runs through `imaging` on a 15 mm sample represented by 1000 points, with derivative option `{'period',7}`. Diffusion and flow are zero; the initial and detected spin states are `Lz` and `L+`, respectively, and the gradient amplitude is 0.015 T/m.

## Sequence and signal processing

The sequence uses zero offset, sweeps [100, 5000] Hz, acquisition sizes [32, 2048], and zero-fill sizes [128, 8192]. Its smoothed saltire chirp has a 20-degree flip angle, 0.015 s duration, 10000 Hz bandwidth, 1000 pulse points, and smoothing factor 20. The coherent steps are set from the reciprocal sweep values, with `delta` equal to one quarter of the first step.

The output FID is reduced to the pure-shift signal by reshaping its first `5000/100` rows. The script applies Gaussian apodisation with parameter 6, calculates a zero-filled 2D FFT of the full FID and a 1D FFT of the reconstructed signal, and plots the 2D magnitude alongside the imaginary part of the 1D spectrum. The source header estimates hours of calculation and notes faster execution on a GPU; GPU enablement is commented out.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
