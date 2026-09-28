# examples/dnp_sol/crosspol_powder_static_1.m

- Signature: `crosspol_powder_static_1()`

## Purpose

Simulates a static-powder (^{15}mathrm{N})–electron cross-polarization contact experiment in the doubly rotating frame and plots the real (^{15}mathrm{N}) (S_x) signal against contact-pulse duration. The source estimates the calculation time as seconds.

## Physical / mathematical content

- The system contains (^{15}mathrm{N}) and an electron at 9.394 T and 298 K. The listed Zeeman scalars are 0 and 2.0023193043622, and the coordinates place the spins 10.05 source-coordinate units apart along (z).
- The simulation uses a 100-step contact sequence. Each step is (10,mumathrm{s}); the electron and nitrogen irradiation-power arrays are both set to (5 × 10^4) for all steps.
- The detected operator is the nitrogen (S_x) state. The source requests an isotropic-equilibrium term and uses the `rep_2ang_6400pts_sph` powder grid.

## Numerical / algorithmic content

Creates the spin system with the full `sphten-liouv` basis (`approximation='none'`) and calls `powder` with the `cp_contact_hard` sequence and NMR mode. The time axis starts at zero and is formed from the cumulative step durations; the plotted signal is the real part of the FID.

## Implementation structure

The function specifies the field, isotopes, Zeeman scalars, coordinates and temperature; constructs and bases the Spinach system; sets the electron and nitrogen irradiation and excitation operators, detected nitrogen state, powder grid and 100-step timing; runs the powder simulation; and plots the resulting nitrogen signal versus contact duration.
