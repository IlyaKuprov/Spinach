# examples/nmr_metabol/molecule_a.m

- Signature: `molecule_a()`

## Purpose

Simulate a ¹H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds.

## Physical / mathematical content

- Imports the spin system from `molecule_a.xml` using `gissmo2spinach('molecule_a.xml',1)`.
- Simulates liquid-state ¹H NMR acquisition with `liquid(spin_system,@acquire,parameters,'nmr')`.

## Numerical / algorithmic content

- Uses the `sphten-liouv` basis formalism with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- Sets both the initial state and detection coil to the ¹H `L+` state; specifies no decoupling.
- Sets offset `3500`, sweep `3000`, `4096` acquisition points, `16536` zero-filled points, ppm axis units, and an inverted axis.
- Applies Gaussian apodisation with parameter `10`, then computes `fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

1. Import the GISSMO dataset and create the Spinach spin system and basis.
2. Set sequence parameters and acquire the FID.
3. Apodise and Fourier-transform the FID.
4. Plot the real spectrum with `plot_1d`.
