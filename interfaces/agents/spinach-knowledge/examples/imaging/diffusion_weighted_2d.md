# examples/imaging/diffusion_weighted_2d.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/diffusion_weighted_2d.m

## Experiment

Constructs a 2D phase-encoded diffusion-weighted image using a spatially varying geometric pattern as the diffusion-coefficient field. The source estimates minutes of runtime and notes a Tesla V100 may shorten it; GPU enablement is commented out, so this function does not request a GPU.

## Spin and image model

The model is a single 1H spin with `sys.magnet=5.9` and zero scalar Zeeman shift. The domain is configured as [0.30 0.25] with a [90 108] grid and third-order periodic derivatives; requested image size is [101 105]. The example loads pattern from ../../etc/phantoms/`pattern.mat`, sets `dxx=dyy=1e-3*pattern`, and sets the off-diagonal tensor components and both in-plane flow fields to zero. Initial magnetisation is `Lz`, detection is `L+`, and both spatial profiles are uniform. No relaxation phantom is provided.

## Encoding and output

The phase-encoded sequence is `phase_enc_2d`. Diffusion-gradient amplitudes are [1e-3 1e-3] T/m; readout and phase-encoding amplitudes are 4.3e-3 and 3.8e-3 T/m. The source sets their durations to 2e-3 and 1e-3, and the echo-time parameter to 1e-2; units for these durations are not annotated in this example. The resulting image is shown next to the loaded coefficient phantom. These are configured simulation parameters, not reported measured image values.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
