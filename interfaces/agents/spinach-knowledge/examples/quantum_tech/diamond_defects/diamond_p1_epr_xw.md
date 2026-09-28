# examples/quantum_tech/diamond_defects/diamond_p1_epr_xw.m

- Signature: `diamond_p1_epr_xw()`

## Purpose

Field-swept powder EPR spectra of a P1 centre in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

The example builds the P1 centre with orientation `111` and `14N`, then calculates electron-spin powder EPR spectra at X and W bands.

## Numerical / algorithmic content

Spinach uses the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`) with powder averaging on `rep_2ang_100pts_sph`. The field sweeps have 1024 points and RSPT order `Inf`; the line width is `1e-4 T`, integration tolerance is `1.0`, and transition-moment tolerance is `0.1`.

## Implementation structure

- Set the P1 orientation and nitrogen isotope, then build the model with `diamond_p1`.
- Set the magnet field to 1 T, construct the Zeeman Hilbert basis, and run Spinach housekeeping.
- Run X-band `fieldsweep` at 9.5 GHz over 0.33–0.35 T and W-band `fieldsweep` at 94 GHz over 3.348–3.36 T.
- Plot both spectra against their returned magnetic-field axes.
