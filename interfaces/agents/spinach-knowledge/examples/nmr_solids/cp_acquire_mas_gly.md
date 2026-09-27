# examples/nmr_solids/cp_acquire_mas_gly.m

- Signature: `cp_acquire_mas_gly()`

## Purpose

Simulates 1H-to-13C cross-polarisation followed by acquisition under MAS in alpha-glycine powder. The source states that the reduced Liouville-space calculation includes correlations up to three spins and takes minutes on a Tesla A100, much longer on a CPU.

## Model and calculation

- Builds the spin system from the glycine log file with `g2spinach`, uses a 9.4 T field, and sets the alpha-glycine isotropic shifts to 176.4, 43.6, 2.6, 3.8, 8.0, 8.0, and 8.0 ppm. Spin temperature is 298 K.
- Uses the `sphten-liouv` basis, approximation `IK-0`, and inter-level 3; enables the greedy option, disables `pt`, and neglects interactions below `2*pi*200` rad/s.
- CP/MAS settings: rotor rate 10,000 Hz, axis `[sqrt(2/3) 0 sqrt(1/3)]`, maximum rank 5, grid `rep_2ang_100pts_sph`, offsets `[2e3 10e3]` Hz, high-power field 83 kHz, CP powers `[60 50]` kHz, and CP duration `50e-5` s. Acquisition uses 50 kHz sweep, 512 points, and 4096-point zero filling.
- Calls `singlerot` with `@cp_acquire_soft`, detects on 13C, then applies exponential apodisation (6), Fourier transforms, and plots the real spectrum.
