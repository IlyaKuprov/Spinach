# examples/nmr_metabol/molecule_c.m

- Signature: `molecule_c()`

## Purpose

Simulates the 1H NMR spectrum of a molecule from the GISSMO database. The source notes a calculation time of seconds.

## Physical / mathematical content

- Imports the spin system and interactions from `molecule_c.xml` using `gissmo2spinach('molecule_c.xml',1)`.
- Simulates a liquid-state 1H NMR acquisition with no decoupling.

## Numerical / algorithmic content

- Uses the `sphten-liouv` basis formalism with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 1.
- Sets both the initial state and detection coil to `L+` for `1H`; uses an offset of 2500, sweep of 5000, 4096 points, and zero filling to 16536 points. The plot axis is in ppm and inverted.
- Applies Gaussian apodisation with parameter 10, then computes `fftshift(fft(fid,parameters.zerofill))` and plots the real spectrum.

## Implementation structure

1. Import the GISSMO dataset and construct the Spinach spin system and basis.
2. Set sequence parameters and acquire the FID with `liquid(spin_system,@acquire,parameters,'nmr')`.
3. Apodise and Fourier-transform the FID, then plot with `plot_1d`.