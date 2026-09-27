# examples/relaxation_theory/dd_csa_xcorr_2.m

- Signature: `dd_csa_xcorr_2()`

## Purpose

DD–CSA cross-correlation example reproducing Fig. 5a from Grace and Kumar ([http://dx.doi.org/10.1006/jmra.1995.1151](http://dx.doi.org/10.1006/jmra.1995.1151)). Calculation time: seconds.

## Imported system and relaxation model

The source reads the vacuum-DFT spin system from `../standard_systems/fdnb.log` using `gparse` and `g2spinach` with the H/`1H` and F/`19F` mappings and arguments `[32.0 270.0]`. It sets the field to 9.4 T, Redfield relaxation with secular retention, Di Bari equilibrium at 298 K, and `tau_c={9.6e-12}`. The proximity cutoff is first assigned 5 Å and then overwritten with 4.0 Å before system creation. The basis is `sphten-liouv` without approximation.

## Simulation and spectra

The acquisition uses `19F`, offset −521, sweep 50, 128 points, and zero-fill 512. The code assumes NMR conditions, forms the Hamiltonian plus `1i*relaxation(spin_system)`, and applies the frequency offset. For each mixing time `[0.1 1.4 1.6 1.8 2.0 2.2 2.4 10]` s, it starts from thermal equilibrium, applies a pi pulse, evolves through mixing, applies a pi/2 pulse, and acquires the detection period. The FID is exponentially apodised with 6 before Fourier transformation; the real spectra are plotted against 19F linear frequency in Hz.
