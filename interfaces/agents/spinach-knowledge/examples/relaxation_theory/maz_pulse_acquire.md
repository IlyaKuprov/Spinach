# examples/relaxation_theory/maz_pulse_acquire.m

- Signature: `maz_pulse_acquire()`

## Purpose

Simulates a proton pulse-acquire spectrum of methylaziridine, illustrating scalar relaxation of the second kind from rapid quadrupolar relaxation of the 14N nucleus. The source estimates minutes of calculation time.

## Physical / mathematical content

The eight-spin model contains seven protons and one 14N nucleus at 11.75 T. Vacuum-DFT shielding tensors, a 14N quadrupole tensor, scalar couplings, and angstrom-scale Cartesian coordinates define the system; isotropic shifts are assigned from experiment. Redfield and SRSK relaxation are combined with zero equilibrium, secular retention, a 200 ps correlation time, and nucleus 4 as the SRSK source.

## Numerical / algorithmic content

The basis is `sphten-liouv` with IK-2 approximation, scalar-coupling connectivity, and proximity level 3. Inter-spin and proximity cutoffs are 2.0 and 4.0, and Krylov propagation is disabled. The proton acquisition starts and detects with `L+`, uses no decoupling, offset 500 Hz, sweep width 1400 Hz, 4096 points, and zero-fills to 16536. Exponential apodisation with parameter 6 precedes the Fourier transform; the plotted spectrum uses the real part and inverted-axis setting.

## Implementation structure

After setting molecular parameters and relaxation, the script builds the spin system and basis, fills the acquisition parameters, and calls `liquid` with `@acquire` in the NMR context. It apodises the FID, Fourier transforms it with `fftshift`, and displays the real spectrum with `plot_1d`.
