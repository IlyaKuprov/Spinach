# examples/nmr_metabol/molecule_b.m

- Signature: `molecule_b()`

## Purpose

Simulate a 1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds.

## Implementation

- Import the GISSMO dataset with `gissmo2spinach('molecule_b.xml',1)`.
- Build the Spinach basis using `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 1.
- Create the spin system and basis with `create` and `basis`.
- Set both the initial state and detection coil to `state(spin_system,'L+','1H')`; use no decoupling.
- Set offset to 3500, sweep to 5000, acquisition points to 4096, zero filling to 16536, axis units to `ppm`, and axis inversion to 1.
- Acquire the liquid-state NMR FID with `liquid(spin_system,@acquire,parameters,'nmr')`, apply Gaussian apodisation with parameter 10, and compute `fftshift(fft(fid,parameters.zerofill))`.
- Plot the real spectrum with `plot_1d`.

Source attribution: ilya.kuprov@weizmann.ac.il