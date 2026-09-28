# examples/nmr_liquids/ct_hsqc_strychnine.m

- Signature: `ct_hsqc_strychnine()`

## Purpose

CT HSQC spectrum of strychnine with natural content of 13C isotope. Calculation time: hours.

## Physical / mathematical content

- Two-dimensional constant-time HSQC of strychnine, simulated for 13C and 1H spins. The spin system is generated with the `13C` and `1H` isotopes, then diluted over 13C isotopomers.
- The FID is apodised with squared-cosine windows; positive and negative components are combined as a States signal before Fourier transformation in both dimensions.

## Numerical / algorithmic content

- Uses a sphten-liouv basis with IK-2 approximation, scalar-coupling connectivity, proximity level 1, and greedy algorithmic settings with `prox_cutoff=4.0`. The calculation is parallelised over 13C isotopomers with `parfor`; no GPU execution is present in this example.
- Sequence settings are `J=140`, sweep `[10000 3000]`, offset `[4000 1000]`, `npoints=[256 256]`, and `zerofill=[512 512]`; the F2 13C channel is decoupled.

## Implementation structure

- Set the strychnine spin system to 5.9 T and construct the selected basis.
- Generate and iterate over 13C isotopomers in parallel, building each basis and simulating CT-HSQC.
- Apodise positive and negative FIDs, form the States signal, Fourier transform both dimensions, and plot the real spectrum.
