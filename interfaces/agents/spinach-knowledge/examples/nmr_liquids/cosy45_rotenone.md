# examples/nmr_liquids/cosy45_rotenone.m

- MATLAB implementation: [examples/nmr_liquids/cosy45_rotenone.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/cosy45_rotenone.m)

- Signature: `cosy45_rotenone()`

## Purpose

Calculates and plots a magnitude-mode 45-degree COSY spectrum for a 22-proton rotenone spin system. The source estimates a calculation time of minutes and cites [DOI 10.1002/jhet.5570250160](https://doi.org/10.1002/jhet.5570250160).

## Spin system and COSY pathway

The model is proton-only at 5.9 T. Its 22 isotropic shifts, in spin order and ppm, are [6.72, 6.40, 4.13, 4.56, 4.89, 6.46, 7.79, 3.79, 2.91, 3.27, 5.19, 4.89, 5.03, 1.72, 1.72, 1.72, 3.72, 3.72, 3.72, 3.76, 3.76, 3.76]. The explicitly assigned scalar couplings, in Hz, are J(3,4)=12.1, J(4,5)=3.1, J(3,5)=1.0, J(3,8)=1.0, J(1,8)=1.0, J(6,7)=8.6, J(5,8)=4.1, J(7,9)=0.7, J(7,10)=0.7, J(9,10)=15.8, J(10,11)=9.8, J(9,11)=8.1, J(13,14)=1.5, J(12,14)=0.9, J(13,15)=1.5, J(12,15)=0.9, J(13,16)=1.5, and J(12,16)=0.9.

The driver uses an IK-2 sphten-Liouville basis with scalar-coupling connectivity, proximity level 1, cutoff 4.0 and greedy selection. It declares three S3 symmetry groups for spins [14-16], [17-19] and [20-22], then runs the COSY sequence through liquid in NMR mode. The COSY implementation starts from proton longitudinal magnetisation, applies a 90-degree first pulse, evolves in F1, selects +1 proton coherence, then applies the requested second-pulse angle, pi/4. It acquires the F2 signal with a proton L+ detection operator. The driver sets no separate receiver phase or additional phase cycle.

## Acquisition and display

The proton offset is 1200 Hz and sweep width is 2000 Hz. The acquired grid is [512, 512] points, zero-filled to [2048, 2048], with axes labelled in ppm. Cosine apodisation precedes a two-dimensional Fourier transform; the plotted quantity is abs(spectrum), so the display is magnitude-mode rather than phase-sensitive.

## Scope

This is a reduced proton-only spin model using the stated couplings, basis approximation and symmetry declarations. The script displays a simulated spectrum and cites a source paper, but contains no experimental spectrum overlay or numerical comparison; it does not itself establish agreement with measurement.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
