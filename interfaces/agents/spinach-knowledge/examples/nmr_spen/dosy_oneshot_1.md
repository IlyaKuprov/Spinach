# examples/nmr_spen/dosy_oneshot_1.m

- Signature: `dosy_oneshot_1()`

## Purpose

Runs a one-shot DOSY imaging sequence for three coupled 1H spins with different longitudinal and transverse relaxation rates. The source estimates minutes on an NVIDIA Tesla A100 and longer on CPU.

## Spin system and relaxation

At 11.7428 T, the shifts are 4.70, 3.50, and 1.50 ppm; couplings are 15 Hz (1–2), 15 Hz (2–3), and 10 Hz (1–3). The T1 values are [0.1952, 0.2100, 0.2500] s and T2 values [0.1602, 0.1802, 0.1902] s. Relaxation uses the diagonal-preserving T1/T2 model with zero equilibrium. The basis is full `sphten-liouv`; perturbation theory is disabled and greedy computation enabled.

## Spatial acquisition

The 15 mm sample is represented by 5000 points with a 7-point periodic derivative stencil. A uniform spatial phantom supplies the relaxation, initial 1H `Lz` state, and receive profile. Acquisition uses 1H, 5 kHz sweep, 1024 points, 32768-point zero-fill, and a 2497.78 Hz offset. The imaging call runs `dosy_oneshot` with gradient amplitude 0.255 T/m, `kappa` 0.2, 1 ms gradient duration, 50 ms diffusion delay, and 0.5 ms gradient-stabilisation delay; the FID is exponentially apodised with parameter 5 before Fourier transformation.
