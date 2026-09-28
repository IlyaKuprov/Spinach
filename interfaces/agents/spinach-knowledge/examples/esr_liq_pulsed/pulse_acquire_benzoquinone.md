# examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m

- Signature: `pulse_acquire_benzoquinone()`

## Purpose

Simulate pulse-acquire FFT ESR of the 2-methoxy-1,4-benzoquinone radical in the liquid state, set to reproduce Figure 1 in http://dx.doi.org/10.1002/mrc.1260280313. A common linewidth is represented by a damping relaxation model. Calculation time: seconds.

## Physical / mathematical content

- The spin system contains one electron and six protons at a magnet induction of 0.33 T. The electron's scalar Zeeman value is 2.004577; its six proton couplings are specified in mT and converted to Hz with `mt2hz`.
- Relaxation uses diagonal damping at a rate of 1e6, with zero equilibrium. The first three protons form an `S3` symmetry group in the basis specification.

## Numerical / algorithmic content

- The simulation uses `liquid` with the `acquire` sequence in ESR mode. The initial state and detection operator are both the electron `L+` state.
- Acquisition specifies a −1e7 Hz offset, a 3e7 Hz sweep, and 1024 points. No apodisation is applied; the FID is Fourier transformed with 4096-point zero filling and `fftshift`.
- The real part of the spectrum is plotted with a GHz lab-frame axis; the parameters request a derivative spectrum and inverted axis.

## Implementation structure

- Define the magnet induction, isotope list, scalar Zeeman value, and electron–proton couplings.
- Configure damping relaxation and a `sphten-liouv` basis without approximation, then create the Spinach spin system and basis.
- Set acquisition and plotting parameters, simulate the FID, apply no apodisation, Fourier transform it, and plot the real spectrum.
