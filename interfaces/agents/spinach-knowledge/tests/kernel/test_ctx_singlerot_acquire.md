# tests/kernel/test_ctx_singlerot_acquire.m

- Signature: `result=test_ctx_singlerot_acquire()`

## Purpose

Exercises `singlerot()` with `acquire()` for a compact anisotropic one-spin MAS calculation, including rotor-space projection.

## Physical / mathematical content

The fixture has one anisotropic `1H` spin at 14.1 T, Zeeman principal values `[-2 -2 4]`, and zero Euler angles. It uses the `sphten-liouv` formalism, no approximation, and projection `+1`.

## Numerical / algorithmic content

The acquisition uses `rho0=coil=L+`, zero offset, 2000 Hz sweep, three points, rotor rate 1000, axis `[1 1 1]`, `max_rank=1`, and the `single_crystal` grid at orientation `[0 0 0]`. The test calls `singlerot(spin_system,@acquire,parameters,`nmr`)`.

## Checks

The FID must contain three finite samples. Its first value must equal the initial coil overlap within absolute and relative tolerance `1e-12`.
