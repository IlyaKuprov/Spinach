# examples/nmr_liquids/gcosy90_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/gcosy90_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/gcosy90_strychnine.m)

## Model and sequence

This wrapper simulates Horne-Morris gradient-selected COSY for strychnine using the 22-spin, proton-only network from `strychnine({'1H'})`; no carbon or nitrogen spins are selected. It sets `sys.magnet=5.9` (field value; no unit is stated here) and calls `liquid(...,@gcosy,parameters,'nmr')`; the gradient-selection pulse program is in `experiments/nmr_liquids/gcosy.m`, not encoded by the wrapper. The selected pathway is `P+N`, which the sequence source describes as the P- and N-type components for echo/anti-echo recombination. The observed channel is `1H`.

## Basis and acquisition

The basis is `sphten-liouv` / `IK-2`, with `scalar_couplings` connectivity and proximity level 1; greedy basis construction uses `prox_cutoff=4.0`. The second pulse angle is `pi/2` rad. Acquisition settings are offset 1200 (unit not specified), sweep 2200 Hz, `npoints=[512 512]`, and `zerofill=[2048 2048]`; displayed axes use ppm. The gradient settings are amplitude 3 Gauss/cm, duration `2e-3` s, stabilisation delay `2e-4` s, and active sample length 1.5 cm. The source estimates calculation time as minutes. No relaxation theory or rates are configured by this example.

## Processing and plot

The wrapper applies squared-cosine apodisation to the positive and negative pathway FIDs, Fourier-transforms F2, combines them as `f1_pos+conj(f1_neg)`, then Fourier-transforms F1. It plots `abs(spectrum)` with positive contours.
