# tests/kernel/test_ctx_floquet_acquire.m

- Signature: `result=test_ctx_floquet_acquire()`

## Purpose

Runs a three-point anisotropic one-spin MAS acquisition through `floquet()` and `acquire()` to exercise the Floquet-context route.

## Physical / mathematical content

The test uses an anisotropic Zeeman tensor with principal values `[-2 -2 4]` and zero Euler angles. Floquet propagation represents the periodic rotor-driven dynamics in the Floquet space.

## Numerical / algorithmic content

The one-spin `1H` system is built at 14.1 T in the `sphten-liouv` formalism with no approximation and projection `+1`. The acquisition uses `rho0=coil=L+`, zero offset, 2000 Hz sweep, three points, rotor rate 1000, axis `[1 1 1]`, `max_rank=1`, and grid `leb_2ang_rank_5`. The test calls `floquet(spin_system,@acquire,parameters,'nmr')`.

## Checks

The test requires three output points, compares the first point with the initial coil overlap to absolute and relative tolerance `1e-12`, and verifies every FID sample is finite.