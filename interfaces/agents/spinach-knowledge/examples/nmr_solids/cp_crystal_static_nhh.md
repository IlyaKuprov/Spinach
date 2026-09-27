# examples/nmr_solids/cp_crystal_static_nhh.m

- Signature: `cp_crystal_static_nhh()`

## Purpose

Simulates static, single-crystal ¹H→¹⁵N cross-polarisation for a ¹⁵N coupled to eight protons. The example retains the full Liouville space; its source notes that this is needed for the non-powder calculation with all spins interacting. The source estimates minutes on a Tesla A100 and longer on a CPU.

## Physical / mathematical content

The nine-spin model specifies isotropic shifts, explicit coordinates for the proton bath around ¹⁵N, and temperature 298 K. The source describes the eight-proton bath as scattered on a 2 Å-radius sphere around ¹⁵N. It samples one crystal orientation, `[pi/3 pi/4 pi/5]`, rather than averaging a powder. The experiment uses 50 kHz spin-lock fields on both channels and detects the ¹⁵N transverse signal during CP.

## Numerical / algorithmic content

The source uses `sphten-liouv` with no basis approximation and calls `crystal` with `@cp_contact_hard`. It requests `aniso_eq` for ¹⁵N and uses 100 time steps of 10 μs. Although a source comment says a GPU is needed, the active setting is `sys.enable={'greedy'}`; `'gpu'` appears only as a commented alternative.

## Implementation structure

- Defines the ¹H₈–¹⁵N spin system, isotropic shifts, coordinates, and temperature.
- Builds the full basis, then creates the CP irradiation operators and ¹⁵N coil state.
- Sets the RF powers, equilibrium requirement, time grid, and single-crystal orientation.
- Runs the crystal CP simulation and plots the real ¹⁵N signal versus time.
