# examples/fitting/nmr_kinetics/glucose_exsy_a.m

- MATLAB implementation: [examples/fitting/nmr_kinetics/glucose_exsy_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/nmr_kinetics/glucose_exsy_a.m)

- Signature: `glucose_exsy_a()`

## Purpose and experiment convention

Fits the source-described 2,2,3,3-tetrafluoroglucose NOESY data by varying a four-pool chemical-exchange model and Redfield relaxation/correlation-time parameters. The function name contains `exsy`, while the source comments call the data NOESY and the simulation explicitly uses `@noesy`; retain that distinction when locating or describing the example. The source estimates minutes but explicitly limits the iteration count, so it supplies no convergence or fit-quality claim.

## Dataset and 36-parameter vector

The objective loads `glucose_expt_a.mat` / `Expression1`, applies `rot90(Expression1,2)/5`, and compares it with the simulated spectrum. The four chemical pools are alpha-inside (spins 1–4), alpha-outside (5–8), beta-inside (9–12), and beta-outside (13–16); all 16 spins are `19F`. The trial-vector roles are:

| Entries | Role |
| --- | --- |
| 1–2 | Forward/reverse entries for the alpha-pool exchange block in `inter.chem.rates` |
| 3–4 | Base-10 exponents used to form the two alternating correlation-time values repeated across pools |
| 5 | Overall spectrum scale (applied with a leading minus sign) |
| 6–21 | Per-spin chemical-shift adjustments to the source's fixed reference values |
| 22–31 | Shared scalar-coupling adjustments for the alpha and beta pool spin groups |
| 32 | Offset in the initial alpha/beta concentration vector passed to `equilibrate` |
| 33–34 | Forward/reverse entries for the beta-pool exchange block |
| 35–36 | Common per-spin `r1_rates` and `r2_rates` values |

No units are assigned here to rate, shift, coupling, concentration-offset, or relaxation entries beyond the operations named in the source. The source's initial vector is `[0.9448, 1.6667, -8.4603, -9.0902, 1.0074, 0.1052, 0.1032, 0.0787, 0.1125, 0.1076, 0.1273, 0.0636, 0.1367, 0.0952, 0.1171, 0.0818, 0.1226, 0.1223, 0.1079, 0.1548, 0.0500, -5.8355, -17.7155, 16.7992, -51.4794, 52.2778, -4.8695, -38.2804, 32.1033, 2.8853, -38.9973, 0.9992, 0.8724, 1.3650, 0.6937, 26.6705]`; these are seed values, not fitted estimates.

## Physical and numerical setup

The field is 9.4; the 16 fluorines are partitioned into those four chemical subsystems. The basis is `sphten-liouv` with no approximation. Relaxation is configured as Redfield plus `t1_t2`, with zero equilibrium and secular retention. The NOESY sequence uses mixing-time value 0.5, offset -49000, sweep [8000 8000], [256 256] points, [1024 512] zero fill, ppm axes, and a chemical-equilibrium `Lz` initial state. The source does not state a unit for the mixing-time value.

`fminsearch` uses the hard cap `MaxIter=10`, `MaxFunEvals=Inf`, iterative display, and `UseParallel=true`. Each objective call simulates the two quadrature FIDs, applies squared-cosine apodisation in both dimensions, performs the 2D transforms, takes the negative scaled real spectrum, and returns the scaled sum of squared matrix residuals. The loaded experimental matrix is transformed before objective evaluation; no plot-only denoising is shown for this dataset.

## Entry-point output

`glucose_exsy_a()` displays the optimiser's final vector but declares no MATLAB output argument. On non-worker evaluations it plots simulated and experimental spectra in adjacent panels using `plot_2d`; the experimental display is labelled positive, and the simulated display both. The script has no figure-save call. The capped search and displayed vector are not evidence of a converged fit.
