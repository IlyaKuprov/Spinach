# examples/nmr_liquids/ecosy_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/ecosy_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ecosy_strychnine.m)

## Model and sequence

This example simulates phase-sensitive E.COSY for strychnine on the 22-spin, proton-only network from `strychnine({'1H'})`; no carbon or nitrogen spins are selected. The wrapper sets `sys.magnet=5.9` (field value; no unit is stated here) and delegates the sequence to `liquid(...,@ecosy,parameters,'nmr')`. The sequence implementation is in `experiments/nmr_liquids/ecosy.m`, not in this wrapper; it cites https://doi.org/10.1021/ja00308a042, https://doi.org/10.1063/1.451421, and https://doi.org/10.1016/0022-2364(87)90102-8. Acquisition and detection use the `1H` channel.

## Basis and acquisition

The basis is `sphten-liouv` / `IK-2`, with `scalar_couplings` connectivity and proximity level 1; greedy basis construction uses `prox_cutoff=4.0`. Settings are offset 1200 (unit not specified), sweep 2200 Hz, `npoints=[512 512]`, and `zerofill=[2048 2048]`; displayed axes use ppm. No relaxation theory or rates are configured by this example. The source estimates calculation time as minutes.

## Processing and plot

The wrapper applies squared-cosine apodisation to the cosine and sine FID components in both dimensions. It Fourier-transforms F2, combines the components as `real(f1_cos)+1i*imag(f1_sin)`, then Fourier-transforms F1. It plots `real(spectrum)` with both signs displayed.
