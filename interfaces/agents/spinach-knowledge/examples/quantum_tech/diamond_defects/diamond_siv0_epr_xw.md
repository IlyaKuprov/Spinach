# examples/quantum_tech/diamond_defects/diamond_siv0_epr_xw.m

- Signature: `diamond_siv0_epr_xw()`

## Purpose

Field-swept powder EPR spectra of SiV0 centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

The example models the SiV0 centre using `29Si` with no `13C` nuclei included in the spin system. It calculates the electron-spin powder EPR spectra at X and W bands.

## Numerical / algorithmic content

The calculation uses the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`) and spherical powder grid `rep_2ang_100pts_sph`.  Each sweep has 2048 field points, RSPT order `Inf`, line width `0.001 T`, integration tolerance `1e-4`, and transition-moment tolerance `0.001`.

## Implementation structure

- Set the SiV0 orientation to `111`, choose `29Si`, and set the number of `13C` nuclei to zero; build the system with `diamond_siv0`.
- Set the magnet field to 1 T, construct the Zeeman Hilbert basis, and run Spinach housekeeping.
- Run X-band `fieldsweep` at 9.5 GHz over 0.1–0.5 T and W-band `fieldsweep` at 94 GHz over 3.2–3.5 T.
- Plot each spectrum against its returned magnetic-field axis.
