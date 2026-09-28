# examples/nmr_spen/idosyzs_test_1.m

- Signature: `idosyzs_test_1()`

## Purpose

Models diffusion attenuation during soft pulses in a simplified Zangger–Sterk pure-shift iDOSY sequence. It simulates signal over a gradient-amplitude series and fits a modified Stejskal–Tanner form, including a fitted gradient shift. The method is described in [the cited JMR paper](https://doi.org/10.1016/j.jmr.2019.02.010). The source estimates seconds on an NVIDIA Tesla A100 and longer on CPU.

## Model and sequence

The model is one 1H spin at 11.7426 T with a 4.6 ppm shift, sample length 15 mm, 4000 spatial points, and a 7-point periodic derivative stencil. The reference diffusion coefficient is 18×10⁻¹⁰ m²/s. The soft pulse shape is `gaussian_1000.pk`, sampled at 100 points, with duration 45 ms and phase π; transmitter offset is 2500 Hz. The gradient duration is 2 ms, diffusion delay 100 ms, and Zangger–Sterk gradient amplitude 0.0053 T/m.

## Simulation and fit

Twenty imaging simulations use gradient amplitudes from 0.01 to 0.40 T/m. The normalised intensities are fitted to `A exp[-D·C·(g-g₀)²]`, where the code's Stejskal–Tanner factor C uses the spin, gradient duration, and diffusion delay. The report gives fitted D (scaled by 10⁻¹⁰) and fitted gradient shift g₀. The simulation supplies no relaxation phantom or relaxation operator, and the initial 1H `Lz` state has white spatial margins.
