# examples/imaging/diffusion_weighted_2d.m

- Signature: `diffusion_weighted_2d()`

## Purpose

Simulates a 2D diffusion-weighted image using an arbitrary geometric pattern as the diffusion-coefficient distribution. Runtime is minutes; a Tesla V100 GPU can make it faster.

## Physical / mathematical content

Models a single 1H spin at 5.9 T with a 0.0 chemical shift on a 0.30 × 0.25 spatial domain. The diffusion tensor has dxx = dyy = 1e-3*pattern and dxy = dyx = 0; flow velocities u and v are zero. The initial state is Lz, and detection uses L+.

## Numerical / algorithmic content

Uses a 90 × 108 spatial grid with third-order periodic derivatives and requests a 101 × 105 image. The phase-encoded 2D sequence uses diffusion gradients [1e-3 1e-3] T/m, a readout gradient of 4.3e-3 T/m for 2e-3 s, a phase-encoding gradient of 3.8e-3 T/m for 1e-3 s, and an echo time of 1e-2 s. The spin offset is 0.0, with no decoupling or relaxation operators.

## Implementation structure

The function diffusion_weighted_2d() creates the Spinach spin system in the sphten-liouv formalism with no basis approximation. It loads pattern from ../../etc/phantoms/pattern.mat, assigns uniform initial-state and coil phantoms, and runs imaging(spin_system,@phase_enc_2d,parameters). It plots the diffusion-weighted image beside the diffusion-coefficient phantom. GPU enablement is present as a commented-out setting.
