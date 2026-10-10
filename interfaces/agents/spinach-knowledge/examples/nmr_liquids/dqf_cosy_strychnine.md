# examples/nmr_liquids/dqf_cosy_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/dqf_cosy_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/dqf_cosy_strychnine.m)

## Model and sequence

This wrapper models strychnine using the 22-spin, proton-only network returned by `strychnine({'1H'})`; no carbon or nitrogen spins are selected. The helper supplies proton shifts and scalar couplings. The wrapper sets `sys.magnet=5.9` (field value; no unit is stated here) and calls `liquid(...,@dqf_cosy,parameters,'nmr')` for a phase-sensitive double-quantum-filtered COSY simulation. The experiment routine, not this wrapper, implements the pulse sequence. Its cited sequence sources are https://doi.org/10.1016/0006-291X(83)91225-1 and https://doi.org/10.1021/ja00388a062. The observed channel is `1H`.

## Basis and acquisition

The basis is `sphten-liouv` / `IK-2`, with `scalar_couplings` connectivity and proximity level 1; the system uses greedy basis construction and `prox_cutoff=4.0`. Acquisition settings are offset 1200 (unit not specified), sweep 2200 Hz, `npoints=[512 512]`, and `zerofill=[2048 2048]`; the displayed axes use ppm. The source estimates calculation time as minutes. No relaxation theory or rates are configured by this example.

## Processing and plot

The wrapper applies cosine apodisation separately to both returned FID components in both dimensions, Fourier-transforms F2, forms the States signal as `real(f1_cos)-1i*real(f1_sin)`, then Fourier-transforms F1. It passes `-real(spectrum)` to `plot_2d` with both signs displayed.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
