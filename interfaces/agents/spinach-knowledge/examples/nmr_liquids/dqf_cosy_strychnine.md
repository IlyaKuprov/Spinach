# examples/nmr_liquids/dqf_cosy_strychnine.m

- Signature: `dqf_cosy_strychnine()`

## Purpose

DQF-COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Two-dimensional 1H DQF-COSY simulation of strychnine. The liquid-state sequence selects double-quantum-filtered scalar-coupling correlations.
- Cosine windows are applied to both cosine and sine FID components; the States signal is formed and Fourier transformed along F2 and F1.

## Numerical / algorithmic content

- Uses the sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1, plus greedy settings with `prox_cutoff=4.0`. Sequence parameters are offset `1200`, sweep `2200`, `npoints=[512 512]`, and `zerofill=[2048 2048]` (1H). The simulation is performed once for the full spin system; it does not use the parallel isotope loop or a GPU.

## Implementation structure

- Build the 1H strychnine spin system at 5.9 T and construct the selected basis.
- Simulate DQF-COSY with the stated 1H acquisition settings, apodise the cosine and sine FIDs, form the States signal, Fourier transform both dimensions, and plot the real spectrum.
