# examples/nmr_liquids/hmbc_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/hmbc_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hmbc_sucrose.m)

- Signature: `hmbc_sucrose()`

## Spin system and method

Builds a sucrose spin system from `../standard_systems/sucrose.log` using vacuum-DFT parameters, mapping H to `1H` and C to `13C`; the `g2spinach` call also receives `[31.8 182.1]`. The wrapper sets `min_j=3.0` and `no_xyz=1`, then replaces isotropic shielding entries at spin indices `[1:19 24:30]` with `[94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]`. The source labels these as experimental isotropic shielding values but does not give their units.

The source comment describes natural 13C content; the wrapper generates 13C isotopomers with `dilute`. It sets `sys.magnet=5.9`, enables `zte` and `greedy`, and uses `prox_cutoff=4.0`. The basis is `sphten-liouv` / `IK-2`, connected by scalar couplings with proximity level 1.

## Acquisition and processing

The wrapper passes `J=140`, `delta_b=60e-3`, `sweep=[6000 2500]`, `offset=[5000 900]`, `npoints=[128 128]`, and `zerofill=[512 512]` to `liquid(...,@hmbc,...,'nmr')`. The dimension order is `{'13C','1H'}`, placing `1H` in the direct/F2 dimension; axis units are explicitly set to ppm. Apart from the axis unit, the source does not annotate units for these numerical parameters.

Each isotopomer is built in an `IK-2` basis, simulated, cosine-apodised in both dimensions, Fourier transformed with the specified zero filling, and accumulated. The wrapper contains no explicit relaxation model or relaxation parameters. Pulse-program internals are delegated to `@hmbc` and are not specified here. It plots `abs(spectrum)` through `plot_2d` with arguments `20,[0.05 1.0 0.05 1.0],2,256,6,'positive'`, after `scale_figure([1.5 2.0])`. The source comment estimates calculation time in seconds.
