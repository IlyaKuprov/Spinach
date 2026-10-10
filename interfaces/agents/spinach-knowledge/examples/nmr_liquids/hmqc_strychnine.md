# examples/nmr_liquids/hmqc_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/hmqc_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmqc_strychnine.m)

- Signature: `hmqc_strychnine()`

## Spin system and method

Calls `strychnine({'13C','1H'})`. The helper's requested-species filter leaves 22 `1H` spins and 21 `13C` sites from its strychnine model; the wrapper then generates natural-abundance 13C isotopomers with `dilute`. The helper attributes isotropic shifts and J-couplings to Berger and Braun except the one-bond C18-H18b coupling, which it cites to http://dx.doi.org/10.1016/j.jmr.2014.02.003; coordinates are attributed to the major conformer in http://dx.doi.org/10.1039/C0CC04114A.

The wrapper sets `sys.magnet=5.9`, enables `zte` and `greedy`, and sets `prox_cutoff=4.0`. Its basis is `sphten-liouv` / `IK-2` with scalar-coupling connectivity and proximity level 1.

## Acquisition and processing

It passes `J=140`, `sweep=[10000 3000]`, `offset=[4000 1000]`, `npoints=[256 256]`, and `zerofill=[512 512]` to `liquid(...,@hmqc,...,'nmr')`. The dimension order is `{'13C','1H'}`, with `1H` decoupling in F1 and `13C` decoupling in F2; the direct/F2 channel is therefore `1H`. Axis units are explicitly ppm; the source does not annotate units for the other numerical parameters.

The wrapper builds and simulates each isotopomer in a `parfor` loop, applies cosine apodisation in both dimensions, sums shifted two-dimensional Fourier transforms, and plots the magnitude spectrum after `scale_figure([1.5 2.0])`, using `plot_2d` arguments `20,[0.05 0.5 0.05 0.5],2,256,6,'positive'`. No explicit relaxation model or relaxation parameters are set in this wrapper. Pulse-program internals are delegated to `@hmqc`. The source comment estimates calculation time in minutes.
