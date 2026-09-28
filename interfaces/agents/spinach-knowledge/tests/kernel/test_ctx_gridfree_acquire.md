# tests/kernel/test_ctx_gridfree_acquire.m

- Signature: `result=test_ctx_gridfree_acquire()`

## Purpose

Exercises `gridfree()` with `acquire()` on a compact anisotropic one-spin MAS model. The test targets the grid-free Fokker–Planck context and its SLE-space projection.

## Physical / mathematical content

The one-spin Zeeman tensor has principal values `[-2 -2 4]` and zero Euler angles. The test uses the `sphten-liouv` formalism, no approximation, and projection `+1`.

## Numerical / algorithmic content

At 14.1 T, the acquisition starts from and detects `L+` on `1H`, with zero offset, 2000 Hz sweep, three points, rotor rate 1000, axis `[1 1 1]`, and `max_rank=2`. It calls `gridfree(spin_system,@acquire,parameters,`nmr`)`.

## Checks

The test verifies the output length is three, the first FID value equals `coil'*rho0` within absolute and relative tolerance `1e-12`, and all samples are finite.
