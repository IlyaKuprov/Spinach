# examples/esr_liq_pulsed/pulse_acquire_chrysene.m

- Signature: `pulse_acquire_chrysene()`

## Purpose

W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in a non-viscous liquid. Simple common line width is used as a relaxation model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2 group direct product. Calculation time: seconds

## Physical / mathematical content

- Reads a vacuum DFT spin system for the chrysene cation radical, including electron and proton spins and their hyperfine couplings, with coordinate information ignored.
- Uses a magnet induction of 3.5 T and diagonal damping relaxation at a rate of 1e6.
- Simulates a liquid-state ESR pulse-acquire signal with electron raising-operator initial state and detection coil.

## Numerical / algorithmic content

- Uses an unrestricted `sphten-liouv` basis and six proton-pair S2 symmetry groups.
- Acquires 1024 points with a 1e8 sweep and −2e7 offset, then Fourier-transforms the FID with zero filling to 4096 points and `fftshift`. Apodisation is set to `none`.
- Plots the real spectrum with a derivative setting, an inverted axis, and `GHz-labframe` axis units.

## Implementation structure

- Suppress coordinate import and load `../standard_systems/chrysene_cation.log` with `gparse` and `g2spinach`.
- Set the magnetic field, relaxation model, basis, and proton-pair symmetry groups; create the Spinach spin system and basis.
- Set sequence parameters and run `liquid(spin_system,@acquire,parameters,'esr')`.
- Apply no apodisation, Fourier-transform the FID, and plot the real spectrum with `plot_1d`.
