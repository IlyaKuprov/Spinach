# examples/fundamentals/nutation_dist_test.m

- Signature: `nutation_dist_test()`

## Purpose

Demonstrates recovery of an RF-field distribution from a nutation curve when the same coil excites and detects the ensemble. The detected signal is reciprocity-weighted by each member’s RF-field amplitude.

## Simulated nutation curve

The system is one on-resonance 1H spin at 14.1 T, in the sphten-liouv formalism with approximation none. The initial state is Lz and the detection components are Lx and Ly. The RF-frequency grid is 2*pi*linspace(25e3,65e3,201) rad/s. Its normalized bimodal Gaussian density has components with weights 0.8 and 0.2, centered at 2*pi*50e3 and 2*pi*38e3 rad/s, with standard deviations 2*pi*3.0e3 and 2*pi*2.5e3 rad/s.

The simulation uses dt=2e-6 s and npts=256. Each ensemble member evolves under H+b1_freq(n)*Lx. Its complex transverse signal is weighted by both its probability mass and RF frequency before accumulation. The curve receives a 1.9-radian phase, is normalized to unit maximum magnitude, and is given reproducible complex Gaussian noise of scale 2e-3 using rng(1).

## Distribution recovery

The script calls nutation_dist(curve,dt,lambda) with lambda=3e2 for second-derivative Tikhonov regularisation, then plots the source and recovered probability densities against nutation frequency in kHz.
