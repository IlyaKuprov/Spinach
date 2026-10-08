# examples/nmr_spen/conv_test.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/conv_test.m)

- Signature: conv_test()

## Purpose

Examines spatial-grid and finite-difference-stencil sensitivity for a pulsed-field-gradient spin echo, then estimates the diffusion coefficient from the simulated attenuation. The source estimates minutes on an NVIDIA Tesla A100 and substantially longer on CPU.

## Spin system and spatial-gradient scan

The spin system is one 1H at field parameter 11.7426, with a single chemical-shift value of 4.6 and no relaxation. The sample length is 0.015 m; velocity is zero. The initial spatial profile is Gaussian with an Lz spin state, and detection uses a spatially uniform L+ coil. Fifty grid sizes span 1000–10000 points, and periodic derivative stencils of 3, 5, and 7 points are compared. The small-gradient duration parameter is 0.002, the diffusion-delay parameter is 0.050, and the reference diffusion coefficient is 18 × 10⁻¹⁰ m²/s. Twenty gradient amplitudes span 0–0.5 T/m.

The imaging call runs st_ideal for each gradient amplitude. Acquisition uses a sweep setting of 5000, 1024 points, zero filling to 32768, a ppm axis, and offset 2500.

## Observable and comparison

For each grid/stencil pair, the simulated intensities are normalised to the zero-gradient signal. The code computes the diffusion estimate from the logarithm of the nonzero-gradient signals divided by the Stejskal–Tanner factors, using the small-gradient duration and the delay corrected by one third of that duration. It plots the estimate against grid size with a reference line at 18 in units of 10⁻¹⁰ m²/s, and also plots elapsed time and the absolute difference between simulated attenuation and the ideal curve on a logarithmic scale.
