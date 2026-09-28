# examples/nmr_solids/cp_matching_3.m

- Signature: `cp_matching_3()`

## Purpose

Maps the ¹H–¹⁵N Hartmann–Hahn response by scanning both spin-lock powers under MAS. The source estimates a calculation time of hours.

## Physical / mathematical content

The model is a single ¹H–¹⁵N pair with specified shifts and coordinates. At a fixed 10 kHz rotor rate, the example scans each channel from 0 to 50 kHz, starts from ¹H transverse magnetisation, and records the final ¹⁵N signal for every power pair.

## Numerical / algorithmic content

The full `sphten-liouv` basis is used with interaction/proximity cutoffs of 5.0/4.0 and `trajlevel` disabled. A 50×50 grid of power pairs is simulated with `singlerot` and `@cp_contact_hard` at MAS axis `[sqrt(2/3) 0 sqrt(1/3)]`, using `parfor` over the inner scan, `max_rank=3`, the `rep_2ang_200pts_oct` grid, and ten 40 μs steps per run. The signal matrix is displayed as an image and updated after each outer-loop row.

## Implementation structure

- Defines the two-spin system, basis, and transverse operators.
- Configures MAS, the initial ¹H state, ¹⁵N coil, time steps, and powder grid.
- Evaluates all 2,500 power pairs and plots the ¹⁵N signal as a two-dimensional image.
