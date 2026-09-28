# examples/nmr_liquids/dqf_cosy_sucrose.m

- Signature: `dqf_cosy_sucrose()`

## Purpose

DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

- Two-dimensional 1H DQF-COSY simulation of sucrose. Magnetic parameters are initialized from a vacuum DFT log; the liquid-state sequence selects double-quantum-filtered scalar-coupling correlations.
- Cosine windows are applied to both cosine and sine FID components; the States signal is formed and Fourier transformed along F2 and F1.

## Numerical / algorithmic content

- Spin-system generation uses `g2spinach` with `min_j=2.0` and `no_xyz=1` for 1H. The field is 5.9 T; the basis is sphten-liouv / IK-2 with scalar-coupling connectivity and proximity level 1, and greedy settings use `prox_cutoff=4.0`. Acquisition settings are offset `800`, sweep `1700`, `npoints=[512 512]`, and `zerofill=[2048 2048]`.

## Implementation structure

- Build the 1H sucrose spin system from the vacuum DFT log, set the 5.9 T field, and construct the selected basis.
- Simulate DQF-COSY with the stated acquisition settings, apodise the cosine and sine FIDs, form the States signal, Fourier transform both dimensions, and plot the real spectrum.
