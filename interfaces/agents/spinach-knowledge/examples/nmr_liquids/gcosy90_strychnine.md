# examples/nmr_liquids/gcosy90_strychnine.m

- Signature: `gcosy90_strychnine()`

## Purpose

Gradient-selected COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Two-dimensional gradient-selected COSY simulation of strychnine using 1H spins and a 90-degree pulse. Scalar-coupling evolution generates the COSY correlations; gradient selection uses the `P+N` pathway.
- The positive and negative echo FIDs are squared-cosine apodised, Fourier transformed along F2, combined as an echo/anti-echo signal, and transformed along F1.

## Numerical / algorithmic content

- Uses the sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1, and greedy settings with `prox_cutoff=4.0`; the field is 5.9 T. Sequence settings are angle `pi/2`, offset `1200`, sweep `2200`, `npoints=[512 512]`, `zerofill=[2048 2048]`, gradient amplitude `3`, duration `2e-3`, stabilization delay `2e-4`, and `s_len=1.5`. The single simulation has no isotopomer-parallel or GPU loop.

## Implementation structure

- Create the 1H strychnine spin system at 5.9 T and construct the selected basis.
- Run gradient-selected COSY with the listed pulse, gradient, and acquisition settings; apply squared-cosine apodisation, form the echo/anti-echo signal, Fourier transform both dimensions, and plot the real spectrum.
