# examples/esr_liq_pulsed/relaxation_bisnitroxide.m

- Signature: `relaxation_bisnitroxide()`

## Purpose

X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical, using explicit time-domain simulation with a Redfield relaxation superoperator. Parameters from https://doi.org/10.1039/C8CP06819D. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains two electrons and two `14N` nuclei, with Zeeman interactions, electron–electron and electron–nitrogen couplings, and a magnetic field of 0.35 T.
- Relaxation uses the `redfield` model with a correlation time of `4e-10` s, zero equilibrium, and lab-frame relaxation terms.

## Numerical / algorithmic content

- The simulation uses the `sphten-liouv` formalism without basis approximation. It acquires 512 points over a `5e8` Hz sweep with zero offset.
- The acquired signal receives no apodisation and is Fourier transformed with zero filling to 1024 points. The real spectrum is plotted using a `GHz-labframe` axis.

## Implementation structure

- Define the spin system, Zeeman interactions, couplings, magnetic field, basis, and relaxation settings.
- Create the Spinach spin system and set the pulse-acquire sequence parameters.
- Run `liquid(spin_system,@acquire,parameters,'esr')`, apply no apodisation, Fourier transform the signal, and plot the real spectrum.