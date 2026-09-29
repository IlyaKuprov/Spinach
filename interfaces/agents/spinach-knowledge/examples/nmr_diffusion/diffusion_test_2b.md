# examples/nmr_diffusion/diffusion_test_2b.m

Source: [examples/nmr_diffusion/diffusion_test_2b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/diffusion_test_2b.m)

## Model

A two-dimensional diffusion-equation example without spin dynamics. It loads `R1` from `phantom_a.mat` and uses `R1(:)` as the initial field; the source does not define the array's physical meaning or units. The spin setup uses ghost isotope `G`, zero magnet field, empty Zeeman and coupling matrices, and the `sphten-liouv` basis with no approximation. Both flow components are zero; the derivative setting is `{'period',7}`.

## Geometry and spatially varying transport

The sample is `[0.02 0.02] m` on a `[108 90]` grid. Both diagonal diffusion components use the same spatial profile, `5e-5*kron(ones(108,1),linspace(0,1,90).^2)`; the off-diagonal components are zero. Thus the diagonal coefficient varies quadratically from zero to `5e-5` along the second grid dimension. The source does not state its unit; metre coordinates and second-based propagation imply m^2/s. This is a spatially varying scalar diffusion field, not a nonzero cross-diffusion tensor. The derivative setting is periodic as specified by `{'period',7}`.

## Propagation and observable

The generator is built with `v2fplanck(spin_system,parameters)` and `inflate`. `evolution` records `200` steps of `5e-4 s` from `R1(:)` (a total modeled interval of `0.1 s`). Each trajectory column is reshaped to `108 x 90` and displayed with `imagesc`. The source labels the calculation time as minutes.
