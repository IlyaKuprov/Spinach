# examples/quantum_tech/diamond_defects/diamond_p1_13c_epr_xw.m

- Signature: `diamond_p1_13c_epr_xw()`

## Purpose

Field-swept powder EPR spectra of a P1 centre in 13C-enriched diamond at X and W bands. Calculation time: minutes.

## Physical / mathematical content

The P1 model uses a `14N` centre with orientation `111` in 13C-enriched diamond. After constructing the full spin system, the example retains the nitrogen and carbon spins whose isotropic hyperfine couplings exceed 8 MHz, then computes electron-spin powder EPR spectra.

## Numerical / algorithmic content

The remaining spin system is treated in the unapproximated Zeeman Hilbert-space basis (`zeeman-hilb`). Powder averaging uses `rep_2ang_100pts_sph`; each field sweep has 1024 points and RSPT order `Inf`. The common line width is `5e-4 T`, with integration and transition-moment tolerances both `0.1`.

## Implementation structure

- Build the P1/13C model with `diamond_p1_13c`, then prune spins to retain 14N and 13C nuclei with isotropic hyperfine coupling above 8 MHz.
- Set the magnet field to 1 T, construct the Zeeman Hilbert basis, and run Spinach housekeeping.
- Simulate X band at 9.5 GHz over 0.31–0.36 T and W band at 94 GHz over 3.33–3.38 T with `fieldsweep`.
- Plot both spectra against their returned magnetic-field axes.
