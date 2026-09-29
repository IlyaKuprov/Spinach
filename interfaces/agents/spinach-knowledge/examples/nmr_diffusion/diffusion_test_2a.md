# examples/nmr_diffusion/diffusion_test_2a.m

Source: [examples/nmr_diffusion/diffusion_test_2a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/diffusion_test_2a.m)

## Model

A two-dimensional diffusion-equation example without spin dynamics. It loads `R1` from `phantom_a.mat`, then uses `R1(:)` as the initial field. The source does not define the physical meaning or units of that loaded array. The spin setup uses ghost isotope `G`, zero magnet field, empty Zeeman and coupling matrices, and the `sphten-liouv` basis with no approximation. Both flow components are zero; the derivative setting is `{'period',7}`.

## Geometry and transport

The rectangular sample is `[0.02 0.02] m` on a `[108 90]` grid. The diffusion tensor is spatially uniform and isotropic: `dxx = dyy = 5e-5`, with `dxy = dyx = 0`. The source does not state the diffusion-coefficient unit; metre coordinates and second-based propagation imply m^2/s.

## Propagation and observable

The generator is built with `v2fplanck(spin_system,parameters)` and `inflate`. `evolution` records `200` steps of `5e-4 s` from `R1(:)` (a total modeled interval of `0.1 s`). Each trajectory column is reshaped to `108 x 90` and displayed with `imagesc`. The source labels the calculation time as minutes.
