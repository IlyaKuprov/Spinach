# tests/kernel/test_liquid_single_spin_fid.m

- Signature: `result=test_liquid_single_spin_fid()`

## Purpose

Tests the free induction decay (FID) of an isolated, zero-offset `1H` spin in the liquid state.

## Physical / mathematical content

With no Hamiltonian-driven precession or relaxation, the detected transverse magnetisation is constant in time.

## Numerical / algorithmic content

- Builds a one-spin system at a magnetic field of `14.1`, with zero scalar Zeeman interaction, using the `sphten-liouv` formalism and no basis approximation.
- Sets both the initial state and detection coil to `L+`, with zero acquisition offset, no decoupling, a sweep of `1000` Hz, `8` points, and `8` zero-fill points.
- Simulates the FID with `liquid(spin_system,@acquire,parameters,'nmr')` and checks that every point equals the first within absolute and relative tolerances of `1e-12`.

## Outputs

- `result` — regression test result with explanatory messages.