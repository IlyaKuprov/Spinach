# examples/fundamentals/high_spin_system_1.m

- Signature: `high_spin_system_1()`

## Purpose

Simulates a pulse-acquire NMR spectrum for a hypothetical scalar coupling to 235U, illustrating the splitting of proton spectral lines.

## Spin system and acquisition

The system is at 14.1 T and contains 1H, 235U, 1H, and 1H spins. Their scalar shifts are −0.5, 0.0, 2.5, and 1.3 ppm; the specified scalar couplings are J12=100 Hz, J34=50 Hz, and J44=0. The basis is sphten-liouv with approximation none.

Acquisition observes 1H only, starting from L+ and detecting with an L+ coil; no decoupling is applied. The offset is 0, sweep width 3500 Hz, and the FID has 1024 points, zero-filled to 4096. The axis is in ppm and inverted.

## Processing

The script calls liquid with acquire and NMR mode, applies exponential apodisation with parameter 6, and computes fftshift(fft(fid,4096)). It plots the real spectrum with plot_1d.
