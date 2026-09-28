# examples/nmr_solids/cp_matching_2.m

- Signature: `cp_matching_2()`

## Purpose

Sweeps proton spin-lock power to test the Hartmann–Hahn matching condition for ¹H–¹⁵N cross-polarisation with a low-power ¹⁵N spin lock. The source describes matching-condition reflections with opposite phase and estimates a calculation time of seconds.

## Physical / mathematical content

The two-spin model uses the specified isotropic shifts and coordinates under a 10 kHz MAS rate. The simulation starts from ¹H transverse magnetisation, scans ¹H power from 0 to 30 kHz while holding ¹⁵N power at 1 kHz, and records the final ¹⁵N signal.

## Numerical / algorithmic content

The full `sphten-liouv` basis is used. With MAS axis `[sqrt(2/3) 0 sqrt(1/3)]`, fifty power points are evaluated in a `parfor` loop with `singlerot` and `@cp_contact_hard`; each run uses `max_rank=3`, the `rep_2ang_200pts_oct` grid, and ten 40 μs steps. The plotted signal is the real part of the final FID point.

## Implementation structure

- Defines the two-spin model and constructs the basis and transverse operators.
- Sets MAS and CP parameters, the ¹H initial state, ¹⁵N coil, time grid, and powder grid.
- Performs the parallel proton-power sweep with the ¹⁵N power fixed at 1 kHz, then plots the ¹⁵N signal.
