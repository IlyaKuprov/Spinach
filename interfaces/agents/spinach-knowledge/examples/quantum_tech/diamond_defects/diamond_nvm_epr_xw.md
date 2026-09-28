# examples/quantum_tech/diamond_defects/diamond_nvm_epr_xw.m

- Signature: `diamond_nvm_epr_xw()`

## Purpose

Field-swept powder EPR spectra of an NV centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

This example models the NV centre using the `14N` isotope and orientation `111`, and computes its electron-spin powder EPR response at X and W microwave frequencies.

## Numerical / algorithmic content

The model uses the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`) and the spherical powder grid `rep_2ang_100pts_sph`. Each field sweep uses 512 points, RSPT order `Inf`, line width `0.001 T`, integration tolerance `0.0001`, and transition-moment tolerance `0.01`.

## Implementation structure

- Set the NV orientation and nitrogen isotope, then construct the ground-state model with `diamond_nvm_gs`.
- Set the magnet field to 1 T, build the Zeeman Hilbert basis, and run Spinach housekeeping.
- Use `fieldsweep` at 9.5 GHz over 0.1–0.5 T (X band) and at 94 GHz over 3.2–3.5 T (W band).
- Plot the two spectra against their returned magnetic-field axes.
