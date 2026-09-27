# examples/esr_liq_pulsed/pulse_acquire_methyl.m

- Signature: `pulse_acquire_methyl()`

## Purpose

X-band pulse-acquire FFT ESR spectrum of methyl radical. A common line width is used as a relaxation model. The example is set to reproduce Figure 4 from the paper by Zhitnikov and Dmitriev: http://dx.doi.org/10.1051/0004-6361:20020268. Calculation time: seconds.

## Physical / mathematical content

- Reads a methyl radical spin system from a vacuum DFT calculation, using supplied hyperfine couplings rather than coordinate information.
- Sets the magnet induction to 0.33 and models relaxation with diagonal damping at a rate of `2.5e7`.
- Simulates an electron-spin pulse-acquire ESR signal with `liquid(spin_system,@acquire,parameters,'esr')`.

## Numerical / algorithmic content

- Uses the `sphten-liouv` basis with no approximation, projection `+1`, and longitudinal `1H` states.
- Sets a sweep of `5e8`, acquires 256 points, and specifies 1024 points for the Fourier transform. The axis units are `GHz-labframe`; the derivative and axis-inversion flags are enabled.
- Applies no apodisation, computes `fftshift(fft(fid,parameters.zerofill))`, and plots the real part of the spectrum.

## Implementation structure

- Ignore coordinate information (HFCs provided).
- Read the spin system (vacuum DFT calculation).
- Set magnet induction, basis, and relaxation parameters.
- Create the Spinach spin system and basis.
- Set the sequence parameters and run the simulation.
- Apply apodisation, perform the Fourier transform, and plot the spectrum.