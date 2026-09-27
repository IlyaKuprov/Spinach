# examples/nmr_liquids/roesy_strychnine.m

- Signature: `roesy_strychnine()`

## Purpose

Simulate and plot a liquid-state ROESY spectrum of strychnine. Stated calculation time: minutes.

## Spin system and settings

- Load strychnine’s `1H` spin system with `strychnine({'1H'})`; set `sys.magnet=5.9`.
- Use the `sphten-liouv` basis with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 3.
- Set Redfield relaxation, zero equilibrium, secular relaxation terms, and correlation time `200e-12`.
- Enable `greedy`, disable `krylov`, and set the proximity cutoff to `4.0`.

## Sequence and processing

- Set mixing time to `0.5`, offset to `1200`, sweeps to `[2500 2500]`, points to `[512 512]`, and zero filling to `[2048 2048]`. Use `1H`, ppm axes, and the `Lz` state for `1H` as the initial state.
- Generate the signal with `liquid(spin_system,@roesy,parameters,'nmr')`. Apply squared-cosine apodisation in both dimensions to the cosine and sine signals.
- Fourier-transform the cosine and sine signals along dimension 1, taking their imaginary and real parts respectively; combine them as `f1_cos-1i*f1_sin`. Fourier-transform the result along dimension 2 and plot the real spectrum with `plot_2d`.