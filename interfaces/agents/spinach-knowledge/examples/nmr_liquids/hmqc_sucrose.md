# examples/nmr_liquids/hmqc_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/hmqc_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmqc_sucrose.m)

- Signature: `hmqc_sucrose()`

## Spin system and method

Builds the sucrose spin system from `../standard_systems/sucrose.log`, mapping H to `1H` and C to `13C` from vacuum-DFT parameters; the `g2spinach` call also receives `[31.8 182.1]`. The wrapper sets `min_j=3.0` and `no_xyz=1`, then replaces isotropic shielding entries at spin indices `[1:19 24:30]` with `[94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]`. The source calls these experimental isotropic shielding values but supplies no units.

The example describes natural 13C content and generates 13C isotopomers with `dilute`. It sets `sys.magnet=5.9`, enables `zte` and `greedy`, and uses `prox_cutoff=4.0`. The basis is `sphten-liouv` / `IK-2`, connected by scalar couplings with proximity level 1.

## Acquisition and processing

It passes `J=140`, `sweep=[10000 3000]`, `offset=[4000 1000]`, `npoints=[256 256]`, and `zerofill=[512 512]` to `liquid(...,@hmqc,...,'nmr')`. The dimension order is `{'13C','1H'}`, with `1H` decoupling in F1 and `13C` decoupling in F2; the direct/F2 channel is therefore `1H`. Axis units are explicitly ppm; the source does not annotate units for the other numerical parameters.

Each isotopomer is simulated in a `parfor` loop, cosine-apodised in both dimensions, Fourier transformed with the specified zero filling, and accumulated. No explicit relaxation model or relaxation parameters appear in this wrapper. Pulse-program internals are delegated to `@hmqc`. It plots `abs(spectrum)` after `scale_figure([1.5 2.0])`, using `plot_2d` arguments `20,[0.05 0.5 0.05 0.5],2,256,6,'positive'`. The source comment estimates calculation time in seconds.
