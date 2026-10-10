# examples/nmr_spen/idosyzs_test_1.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/idosyzs_test_1.m)

- Signature: idosyzs_test_1()

## Purpose

Models diffusion attenuation during soft pulses in a simplified Zangger–Sterk pure-shift iDOSY sequence. The sequence and modified Stejskal–Tanner fit are described in [the cited Journal of Magnetic Resonance paper](https://doi.org/10.1016/j.jmr.2019.02.010). The source estimates seconds on an NVIDIA Tesla A100 and substantially longer on CPU.

## Spin system and spatial selection

The model is one 1H at field parameter 11.7426, with shift value 4.6 and a reference diffusion coefficient of 18 × 10⁻¹⁰ m²/s. The sample length is 0.015 m, represented by 4000 points with a 7-point periodic derivative stencil. Relaxation phantoms and operators are empty. The initial Lz state occupies the central 2000 spatial points with 1000-point zero margins on both sides; detection uses a uniform L+ phantom.

## Sequence, encoding, and fit

The soft-pulse shape is gaussian_1000.pk, sampled at 100 points with duration parameter 0.045; the inversion-pulse phase is pi. The sequence sets transmitter offset 2500, small-gradient-duration parameter 0.002, diffusion-delay parameter 0.1, and Zangger–Sterk selection-gradient amplitude 0.0053. It runs idosyzs for 20 diffusion-gradient amplitudes spanning 0.01–0.40 T/m.

The simulated intensities are normalised to the first point and fitted to an exponential attenuation with an amplitude, diffusion coefficient, and fitted gradient shift. The fit factor uses the spin, gradient duration, and delay corrected by one third of the gradient duration, with the source’s 10⁻¹⁰ scaling; the reported diffusion coefficient is the fitted coefficient multiplied by 10⁻¹⁰. The page plots simulated points and the fitted curve and reports the fitted diffusion coefficient and gradient shift.
