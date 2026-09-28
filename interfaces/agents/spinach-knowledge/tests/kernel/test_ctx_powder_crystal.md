# tests/kernel/test_ctx_powder_crystal.m

- Signature: `result=test_ctx_powder_crystal()`

## Purpose

Compares `powder()` and `crystal()` for a single static orientation represented by the `single_crystal` grid.

## Physical / mathematical content

The fixture is one anisotropic `1H` spin at 14.1 T, with Zeeman principal values `[-2 -2 4]` and zero Euler angles. The selected orientation is `[0 0 0]`.

## Numerical / algorithmic content

The test uses `sphten-liouv`, no approximation, projection `+1`, zero offset, 2000 Hz sweep, four points, and `serial=true`; `rho0` and `coil` are both `L+`. It passes the same parameters and `single_crystal` orientation to both contexts.

## Check

The two FIDs must agree to absolute and relative tolerance `1e-12`. The test identifies the single-point powder grid as unit weight at the zero Euler orientation.
