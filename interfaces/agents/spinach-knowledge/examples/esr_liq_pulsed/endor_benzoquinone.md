# examples/esr_liq_pulsed/endor_benzoquinone.m

- Signature: `endor_benzoquinone()`

## Purpose

Simulates liquid-state continuous-wave ENDOR of the 2-methoxy-1,4-benzoquinone radical, reproducing Figure 2 of http://dx.doi.org/10.1002/mrc.1260280313. Calculation time: seconds.

## Physical / mathematical content

The spin system contains one electron and six protons with scalar electron–proton hyperfine couplings. The basis imposes `S3` symmetry on the three equivalent proton spins; the detected ENDOR channel is proton-based.

## Numerical / algorithmic content

The CW ENDOR calculation uses a 50 MHz sweep, 1024 acquired points, and zero filling to 4096. It subtracts the mean, applies a Kaiser apodisation (parameter 20), Fourier transforms, and plots the negative spectrum magnitude.

## Implementation structure

It creates a full sphten-liouv basis, calls `liquid` with `@endor_cw`, then performs the stated signal processing and plots the result using `plot_1d`.
