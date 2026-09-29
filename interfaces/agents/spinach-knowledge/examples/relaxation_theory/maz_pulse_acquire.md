# examples/relaxation_theory/maz_pulse_acquire.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/maz_pulse_acquire.m) · Signature: `maz_pulse_acquire()`

## Purpose and physical model

This methylaziridine pulse-acquire calculation illustrates scalar relaxation of the second kind from fast quadrupolar relaxation of the `14N` nucleus. It models seven `1H` spins and one `14N` spin at magnet setting 11.75. The source specifies vacuum-DFT shielding tensors, the nitrogen quadrupole tensor, scalar couplings, and molecular coordinates; the coordinate comment states Angstrom. The nitrogen quadrupole matrix is [-1.2932, 0.6251, 1.8700; 0.6251, 1.7170, -2.3127; 1.8700, -2.3127, -0.4238] multiplied by 1e6; the source does not label its units. Scalar-coupling entries from nitrogen to protons include 4.4, 5.2, and 44.8, with units likewise unstated. Its isotropic shielding components are assigned from experiment as 1.3, 1.7, 1.9, 0.0, 0.1, 1.2, 1.2, and 1.2. Those shifts are calculation inputs, not a measured spectrum provided as an output reference. Estimated calculation time is minutes.

## Relaxation and acquisition

The source requests Redfield relaxation together with SRSK, sets the equilibrium state to zero, selects nitrogen spin 4 as the SRSK source, and keeps the relaxation superoperator secular. The correlation-time input is `200e-12` s. The basis is `sphten-liouv` with `IK-2` approximation, scalar-coupling connectivity and proximity level 3; inter-spin and proximity cutoffs are 2.0 and 4.0, and Krylov propagation is disabled.

The source runs `liquid(spin_system,@acquire,parameters,'nmr')`. The `@acquire` simulation observes `1H`, initialises and detects with proton `L+`, and has no decoupling spins configured. It uses offset 500 Hz, sweep 1400, 4096 acquired points, zero-filling to 16536, ppm axis units, and an inverted axis. The FID is exponentially apodised with the source's value 6, Fourier transformed, and the real part of the calculated spectrum is plotted with `plot_1d`. The resulting line shape and intensities are simulated, not a digitised or measured pulse-acquire spectrum.
