# examples/nmr_spen/psyche_rotenone.m

- MATLAB implementation: [examples/nmr_spen/psyche_rotenone.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/psyche_rotenone.m)

- Signature: `psyche_rotenone()`
- Source: [`examples/nmr_spen/psyche_rotenone.m`](../../../../../examples/nmr_spen/psyche_rotenone.m)

## Model and spatial encoding

The example constructs a 22-site `1H` rotenone system with a field of 11.7 T, explicitly assigned chemical-shift entries (offset by 5.0 in the source), and a sparse network of scalar-coupling entries. It uses the spherical-tensor Liouville formalism, IK-2 basis approximation, and scalar-coupling connectivity. The PSYCHE sequence is run through `imaging` on a 15 mm sample represented by 1000 points; the derivative option is `{'period',7}`. Diffusion and flow are set to zero. The initial state is `Lz`, detection is `L+`, and the gradient amplitude is 0.015 T/m.

## Sequence and signal processing

The sequence uses zero offset, sweeps [100, 5000] Hz, acquisition sizes [32, 2048], and zero-fill sizes [128, 8192]. Its smoothed saltire chirp has a 20-degree flip angle, 0.015 s duration, 10000 Hz bandwidth, 1000 pulse points, and smoothing factor 20. The source sets the coherent steps to the reciprocal sweep values and `delta` to one quarter of the first step.

After imaging, the script reshapes the first `5000/100` rows of the FID into a pure-shift FID, applies Gaussian apodisation with parameter 6, and computes a zero-filled 2D FFT of the full data and a 1D FFT of the reconstructed signal. It plots the 2D magnitude and the imaginary part of the 1D spectrum. The source header estimates hours of calculation and says a GPU is faster; GPU enablement is commented out in this example.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
