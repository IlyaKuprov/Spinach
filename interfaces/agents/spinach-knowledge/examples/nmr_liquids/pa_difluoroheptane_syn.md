# examples/nmr_liquids/pa_difluoroheptane_syn.m

- Signature: `pa_difluoroheptane_syn()`

## Purpose

Pulse-acquire a ¹H NMR spectrum of syn-3,5-difluoroheptane. The basis is manually specified by merging Lie algebras of user-selected structural fragments, followed by symmetry factorisation and conservation-law screening. Source paper: [https://doi.org/doi/10.1021/acs.joc.4c00670](https://doi.org/doi/10.1021/acs.joc.4c00670). The source comment estimates minutes of computation and notes it is faster with a GPU.

## Physical / mathematical content

- Models liquid-state NMR with explicit isotopes, chemical shifts, and scalar couplings for syn-3,5-difluoroheptane.
- Uses a fragment-based Liouville basis, then applies symmetry factorisation and conservation-law screening.

## Numerical / algorithmic content

- Constructs the specified basis and simulates a liquid-state acquisition to obtain an FID.
- Applies exponential apodisation, zero-filled Fourier transformation, and plotting of the real spectrum on an inverted ppm axis.

## Implementation structure

- Defines the spin-system parameters and manually assembles the basis; runs the acquisition, apodises and Fourier-transforms the FID, then plots the spectrum.
