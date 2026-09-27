# examples/esr_liq_pulsed/pulse_acquire_phenyl.m

- Signature: `pulse_acquire_phenyl()`

## Purpose

W-band pulse-acquire FFT ESR spectrum of phenyl radical. Simple fixed line width is used as a relaxation model. Calculation time: seconds

## Physical / mathematical content

- Reads phenyl spin-system properties from a vacuum DFT calculation in `../standard_systems/phenyl.log`; coordinate information is ignored because hyperfine couplings are provided.
- Sets the magnet induction to 3.5 and uses diagonal damping relaxation with a rate of `1e7` and zero equilibrium.
- Acquires an electron-spin signal with `L+` as both the initial state and detection state.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism without basis approximation, with `1H` longitudinal states and projection `+1`.
- Runs `liquid(spin_system,@acquire,parameters,'esr')` with a sweep of `2e8`, 512 points, and no decoupling.
- Applies no apodisation, then computes `fftshift(fft(fid,1024))` and plots the real spectrum. The plot parameters request `GHz-labframe` axis units, a derivative, and an inverted axis.

## Implementation structure

- Set `options.no_xyz=1` and read the spin system with `gparse` and `g2spinach`.
- Configure the magnet, basis, relaxation, and spin system.
- Set acquisition parameters and simulate the ESR free-induction decay.
- Apply no apodisation, Fourier-transform the result, and plot the real spectrum.
