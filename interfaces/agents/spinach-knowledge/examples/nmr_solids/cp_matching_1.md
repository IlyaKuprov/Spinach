# examples/nmr_solids/cp_matching_1.m

- Signature: `cp_matching_1()`

## Purpose

Sweeps the proton spin-lock power to examine the Hartmann–Hahn matching condition for ¹H–¹⁵N cross-polarisation under MAS. The source estimates a calculation time of seconds.

## Physical / mathematical content

The model contains one ¹H and one ¹⁵N with the specified isotropic shifts and coordinates. At a 10 kHz rotor rate, the simulation starts from ¹H transverse magnetisation, applies a fixed 50 kHz ¹⁵N spin-lock field, and records the final ¹⁵N signal while the ¹H power varies from 20 to 80 kHz.

## Numerical / algorithmic content

The full `sphten-liouv` basis is used. With MAS axis `[sqrt(2/3) 0 sqrt(1/3)]`, each of 120 power values is simulated in a `parfor` loop with `singlerot` and `@cp_contact_hard`, using the `rep_2ang_200pts_oct` grid, `max_rank=3`, and ten 40 μs time steps. The plotted quantity is the real part of the final FID point.

## Implementation structure

- Defines the two-spin system and builds its basis and transverse operators.
- Constructs the MAS and CP parameters, initial ¹H state, ¹⁵N coil, time grid, and powder grid.
- Runs the parallel single-power sweep and plots ¹⁵N signal against ¹H spin-lock power.
