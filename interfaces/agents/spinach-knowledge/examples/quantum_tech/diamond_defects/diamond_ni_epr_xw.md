# examples/quantum_tech/diamond_defects/diamond_ni_epr_xw.m

- Signature: `diamond_ni_epr_xw()`

## Purpose

Field-swept powder EPR spectra of Ni defects in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

The example builds the diamond Ni NE1 centre with orientation `111` and simulates its electron-spin EPR response. It computes powder field sweeps at X and W microwave bands and plots the two spectra side by side.

## Numerical / algorithmic content

Spinach uses the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`). The field-sweep calculation uses the spherical orientation grid `rep_2ang_100pts_sph`, 2048 field points, and an RSPT order of `Inf`. The common line width is `1e-4 T`; integration tolerance is `1e-4` and transition-moment tolerance is `0.1`.

## Implementation structure

- Set the NE1 centre and `111` orientation, then build the spin system with `diamond_ni`.
- Set the magnet field to 1 T, construct the Zeeman Hilbert basis, and run Spinach housekeeping.
- Simulate X band at 9.5 GHz over 0.30–0.36 T and W band at 94 GHz over 3.1–3.4 T using `fieldsweep`.
- Plot each spectrum against its returned magnetic-field axis.
