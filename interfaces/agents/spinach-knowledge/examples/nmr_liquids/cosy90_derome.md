# examples/nmr_liquids/cosy90_derome.m

- MATLAB implementation: [examples/nmr_liquids/cosy90_derome.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/cosy90_derome.m)

- Signature: `cosy90_derome()`

## Purpose

A three-proton 90-degree COSY simulation identified in the source as Figure 8.26 of Andrew Derome's *Modern NMR Techniques for Chemistry Research*. The driver estimates a calculation time of seconds.

## Spin system and COSY pathway

The system contains three 1H spins at 16.1 T. Their isotropic shifts are 3.70, 3.92 and 4.50 ppm; the explicitly assigned scalar couplings are J(1,2)=10 Hz, J(2,3)=12 Hz and J(1,3)=4 Hz. A zero diagonal coupling is also assigned to spin 3. The driver uses the full sphten-Liouville basis (approximation 'none'), with greedy selection enabled and proximity cutoff 4.0.

The driver runs the COSY pulse program through liquid in NMR mode. It begins with proton longitudinal magnetisation and applies a 90-degree first pulse. It evolves the spin state over F1, selects +1 proton coherence, applies a second 90-degree pulse (angle pi/2), then records F2 using the proton L+ detection operator. No separate receiver phase or additional phase cycle is set in the example driver.

## Acquisition and display

The proton offset is 2800 Hz and sweep width is 700 Hz. The acquired grid is [1024, 1024] points and the data are zero-filled to [2048, 2048]; both axes are labelled in ppm. Square-cosine apodisation is applied in both dimensions before the two-dimensional Fourier transform. The plot uses the real part of the transformed signal.

## Scope

The calculation is an idealised three-spin demonstration tied to the cited book figure, not a molecularly complete spectrum or an experimental data set. The page describes the spin parameters, pulse-program pathway and processing specified by the example; it does not claim a measured-spectrum match.
