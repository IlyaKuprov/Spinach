# examples/quantum_tech/diamond_defects/diamond_ninter_epr_xw.m

- Signature: `diamond_ninter_epr_xw()`

## Purpose

Field-swept powder EPR spectra of a nitrogen interstitial defect in diamond at X and W bands. Calculation time: seconds.

## Physical / mathematical content

The model is the WAR9 nitrogen-interstitial centre in diamond, with `15N` and orientation `111`. The example calculates electron-spin powder EPR field sweeps at X and W bands.

## Numerical / algorithmic content

The spin system is evaluated in the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`). Powder averaging uses `rep_2ang_100pts_sph`; each spectrum has 1024 field points and RSPT order `Inf`. The line width is `1e-5 T`, with integration and transition-moment tolerances both `0.1`.

## Implementation structure

- Set the WAR9 centre, `15N`, and orientation `111`; build the model using `diamond_n_inter`.
- Set the magnet field to 1 T, form the Zeeman Hilbert basis, and run Spinach housekeeping.
- Run X-band `fieldsweep` at 9.755 GHz over 0.347–0.349 T and W-band `fieldsweep` at 94 GHz over 3.351–3.355 T.
- Plot both spectra against their returned field axes.
