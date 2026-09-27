# examples/nmr_solids/cp_crystal_static_nh.m

- Signature: `cp_crystal_static_nh()`

## Purpose

Simulates a static, single-crystal ¹H–¹⁵N cross-polarisation experiment in the doubly rotating frame. The source estimates a calculation time of seconds.

## Physical / mathematical content

The model is a two-spin ¹⁵N–¹H pair with zero isotropic shifts, a 1.05 Å internuclear separation, and temperature 298 K. One specified crystal orientation, `[pi/3 pi/4 pi/5]`, is used; this is not a powder average. The simulated observable is the ¹⁵N transverse signal during CP.

## Numerical / algorithmic content

The source uses the full `sphten-liouv` basis without approximation and calls `crystal` with `@cp_contact_hard`. It requests `aniso_eq` for the ¹⁵N spin, applies 50 kHz spin-lock fields on both spins, and propagates 100 time steps of 10 μs.

## Implementation structure

- Defines the two isotopes, isotropic shifts, coordinates, and temperature.
- Builds the full basis and spin system.
- Sets the two-channel RF operators, ¹⁵N coil state, initial equilibrium requirement, time steps, and crystal orientation.
- Runs the single-crystal CP simulation and plots the real ¹⁵N signal versus time.
