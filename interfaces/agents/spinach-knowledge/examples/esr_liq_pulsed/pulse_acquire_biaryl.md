# examples/esr_liq_pulsed/pulse_acquire_biaryl.m

- Signature: `pulse_acquire_biaryl()`

## Purpose

A time-domain pulse-acquire version of the EasySpin biaryl test file, with acknowledgements to Stefan Stoll. The Spinach simulation uses explicit time propagation in Liouville space and symmetry treatment with the S2xS2xS2xS2xS2xS2 direct-product group. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains one electron, two `14N` nuclei, and ten `1H` nuclei in a 0.33 T magnetic field.
- The source specifies an electron Zeeman scalar of 2.00316, electron–nucleus scalar couplings, and diagonal damping relaxation at a rate of `5e5`.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism without approximation and defines six `S2` symmetry groups over paired nuclei.
- The ESR acquisition uses electron `L+` states for the initial state and detection coil, a `3e8` sweep, and 4096 points. The FID is passed through apodisation set to `none`, then Fourier-transformed with zero filling to 16384 points.

## Implementation structure

- Define the magnet, isotopes, basis, Zeeman interaction, scalar couplings, and relaxation settings.
- Create the spin system and basis, then run `liquid(spin_system,@acquire,parameters,'esr')`.
- Compute `fftshift(fft(fid,parameters.zerofill))` and plot the real spectrum with `plot_1d`.
