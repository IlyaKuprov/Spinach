# examples/nmr_solids/cp_contact_mas_nh_gridfree.m

- Signature: `cp_contact_mas_nh_gridfree()`

## Purpose

Simulates ¹H→¹⁵N cross-polarisation in the doubly rotating frame for a single proton–nitrogen pair. The grid-free Fokker–Planck calculation starts from thermal equilibrium and averages a spinning powder. The source estimates minutes on a Tesla A100 GPU and substantially longer on a CPU.

## Physical / mathematical content

The two-spin system has zero isotropic Zeeman shifts, a 1.05 Å internuclear separation, and temperature 298 K. The experiment applies spin-lock fields to ¹H and ¹⁵N during magic-angle spinning; the rotor axis is `[sqrt(2/3) 0 sqrt(1/3)]`. The detected observable is the ¹⁵N transverse magnetisation.

## Numerical / algorithmic content

The source uses the full `sphten-liouv` basis (`bas.approximation='none'`) and calls `gridfree` with `@cp_contact_hard`. It requests `iso_eq`, uses 100 time steps of 10 μs, and sets `max_rank=42`. The spin-lock powers are 50 kHz on ¹H and 40 kHz on ¹⁵N; the rotor rate is 10 kHz. A source comment says a GPU is needed and shows `sys.enable={'gpu'}` as a commented-out line, so the script does not explicitly enable that option.

## Implementation structure

- Defines the ¹⁵N–¹H pair, isotropic shifts, internuclear coordinates, and temperature.
- Builds the Spinach system and the transverse operators used for irradiation and detection.
- Sets MAS axis, rank, RF powers, equilibrium requirement, and time grid.
- Runs the grid-free CP simulation and plots the real ¹⁵N signal versus accumulated time.
