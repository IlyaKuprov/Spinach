# examples/relaxation_theory/trosy_double.m

- Signature: `trosy_double()`

## Purpose

Hari Arthanari's Double TROSY effect. Calculation time: seconds.

## Physical / mathematical content

- Uses a four-spin system comprising two `1H`, one `19F`, and one `13C` spin, selected from a parsed DFT calculation with their coordinates, chemical-shift tensors, and scalar J-couplings.
- Uses Redfield relaxation with the secular terms retained, zero equilibrium, and a correlation time of `20e-9` s. The magnetic field is `14.1` T.
- Compares the calculated `13C` spectra before and after clearing the coordinates of both protons for a run labeled in the source as having no proton DD.

## Numerical / algorithmic content

- Uses the `sphten-liouv` basis formalism with no basis approximation.
- Initializes and detects `13C` `L+` coherence, with no decoupling. Acquisition uses an offset of `26800`, a sweep of `500`, and `2048` points; the chemical-shift axis is in ppm and inverted.
- Applies Gaussian apodisation with parameter `6`, then computes `fftshift(fft(fid,16384))`. Plots the real spectra in two panels.

## Implementation structure

- Parses `../standard_systems/4_fluoro_phe.out` using `gparse` and `g2spinach`, then selects DFT spin indices `[10 20 19 8]`.
- Creates the spin system and basis, runs `liquid` acquisition, and plots the full spectrum. It then clears the two proton coordinates, rebuilds the spin system, repeats acquisition and processing, and plots the comparison spectrum.
