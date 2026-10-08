# examples/nmr_liquids/dqf_cosy_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/dqf_cosy_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/dqf_cosy_sucrose.m)

## Model and sequence

This wrapper builds the sucrose spin system from the vacuum-DFT file `../standard_systems/sucrose.log` using `gparse` and `g2spinach`, mapping hydrogen atoms to `1H`. The log identifies sucrose as C12H22O11, so the selected network is its 22 proton spins; carbon and oxygen are not selected as spins. The parser options are `min_j=2.0` and `no_xyz=1` (the wrapper does not annotate units for these values). It sets `sys.magnet=5.9` (field value; no unit is stated here) and calls `liquid(...,@dqf_cosy,parameters,'nmr')`. The DQF-COSY pulse program is implemented outside the wrapper in `experiments/nmr_liquids/dqf_cosy.m`, which cites https://doi.org/10.1016/0006-291X(83)91225-1 and https://doi.org/10.1021/ja00388a062. The observed channel is `1H`.

## Basis and acquisition

The basis is `sphten-liouv` / `IK-2`, with `scalar_couplings` connectivity and proximity level 1; greedy basis construction uses `prox_cutoff=4.0`. Acquisition settings are offset 800 (unit not specified), sweep 1700 Hz, `npoints=[512 512]`, and `zerofill=[2048 2048]`; displayed axes use ppm. The source estimates calculation time as minutes. No relaxation theory or rates are configured by this example.

## Processing and plot

Cosine apodisation is applied to both cosine and sine FID components in both dimensions. After the F2 transform, the States signal is formed as `real(f1_cos)-1i*real(f1_sin)` and Fourier-transformed along F1. The plot uses `-real(spectrum)` and displays both signs.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
