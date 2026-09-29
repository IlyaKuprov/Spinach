# examples/nmr_liquids/hmqc_cyprinol.m

- MATLAB implementation: [examples/nmr_liquids/hmqc_cyprinol.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmqc_cyprinol.m)

- Signature: `hmqc_cyprinol()`

## Spin system and method

Loads the helper-defined cyprinol spin system with `cyprinol()`. That helper declares 42 `1H` spins and 27 `13C` sites; the wrapper describes natural 13C abundance and enumerates 13C isotopomers with `dilute`. The helper comments attribute isotropic shifts and J-couplings to http://dx.doi.org/10.1002/mrc.4782, with unspecified values estimated rather than reported in that source.

The wrapper sets `sys.magnet=11.7`, enables `greedy`, and sets proximity and interaction cutoffs to `4.0` and `5.0`. It uses a `sphten-liouv` / `IK-1` basis with interaction level 3, proximity level 1, and scalar-coupling connectivity.

## Acquisition and processing

It passes `J=150`, `sweep=[12000 2500]`, `offset=[5000 1250]`, `npoints=[128 128]`, and `zerofill=[512 512]` to `liquid(...,@hmqc,...,'nmr')`. The dimension order is `{'13C','1H'}`, with `1H` decoupling in F1 and `13C` decoupling in F2; the direct/F2 channel is therefore `1H`. Axis units are explicitly ppm. The wrapper does not annotate units for the other numerical parameters.

Each isotopomer is simulated in a `parfor` loop, cosine-apodised in both dimensions, Fourier transformed and summed. No explicit relaxation model or relaxation parameters appear in this wrapper; pulse-program internals are delegated to `@hmqc`. The absolute spectrum is plotted after `scale_figure([1.5 2.0])` using `plot_2d` arguments `20,[0.05 0.5 0.05 0.5],2,256,6,'positive'`; the source comment estimates calculation time in seconds.
