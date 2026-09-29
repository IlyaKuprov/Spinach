# examples/nmr_diffusion/diffusion_test_1.m

Source: [examples/nmr_diffusion/diffusion_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/diffusion_test_1.m)

## Model

A one-dimensional diffusion-equation example without spin dynamics. It uses a ghost isotope `G`, zero magnet field, empty Zeeman and coupling matrices, and the `sphten-liouv` basis with no approximation. Flow is zero. The derivative setting is `{'period',7}`.

## Geometry and initial profile

The sample length is `0.02 m`, represented by `100` points. The diffusion parameter is `5e-5` (the source does not state its unit; with position in metres and time in seconds, its dimensional unit is m^2/s). The dimensionless initial profile is `exp(-0.125*((1:100)-20).^2)`, a Gaussian-shaped concentration profile centred at grid index 20.

## Propagation and observable

The transport generator is formed with `v2fplanck(spin_system,parameters)` and expanded with `inflate`. `evolution` records `90` trajectory steps at `5e-4 s` each. The script plots concentration against the physical coordinate from `-0.01 m` to `0.01 m`, with the display range fixed to 0–1. The source labels the calculation time as seconds.
