# examples/fitting/nmr_kinetics/glucose_exsy_b.m

- MATLAB implementation: [examples/fitting/nmr_kinetics/glucose_exsy_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/nmr_kinetics/glucose_exsy_b.m)

- Signature: `glucose_exsy_b()`

## Purpose and experiment convention

Fits the source-described 3,3-difluoroglucose NOESY data with exchange and rotational-correlation-time parameters in Redfield theory. As in the companion example, the function name says `exsy`, but the source describes NOESY and calls `@noesy`. Its source comment estimates hours and notes the iteration cap; it gives no measured fit result or convergence claim.

## Dataset and 20-parameter vector

The objective loads `glucose_expt_b.mat` / `spec` and transposes it with `atranspose` before comparing to the simulated matrix. Four chemical pools contain two `19F` spins each: alpha-inside (1–2), alpha-outside (3–4), beta-inside (5–6), and beta-outside (7–8). The fit vector controls:

| Entries | Role |
| --- | --- |
| 1–2 | Forward/reverse entries for the alpha-pool exchange block |
| 3, 14 | Initial alpha and beta concentrations passed to `equilibrate` |
| 4–11 | Per-spin chemical-shift adjustments to fixed reference values |
| 12–13 | Shared scalar-coupling adjustments for alpha and beta spin pairs |
| 15–16 | Forward/reverse entries for the beta-pool exchange block |
| 17–18 | Base-10 exponents used to form the two alternating correlation-time values repeated across pools |
| 19–20 | Common per-spin `r1_rates` and `r2_rates` values |

The initial guess is `[0.2331, 0.1114, 6.8446, 0.1210, 0.1423, 0.1205, 0.1697, 0.0755, 0.1158, 0.0936, 0.1256, 26.9142, 18.7751, 25.5476, 0.7948, 0.4331, -9.0277, -9.2828, 1.1924, 32.9179]`; it is a starting vector, not a reported estimate. The source does not state units for these fitted quantities.

## Physical and numerical setup

The field is 9.3933; the basis is `sphten-liouv` without approximation. Redfield plus `t1_t2` relaxation uses zero equilibrium and secular retention. The NOESY sequence uses mixing-time value 0.5, offset -46681, sweep [8650 8650], [512 1024] points, [1024 1024] zero fill, ppm axes, and an equilibrium-chemical `Lz` initial state. The source does not state a unit for the mixing-time value. It simulates the quadrature FIDs with `@noesy`, applies squared-cosine apodisation in both dimensions, and Fourier transforms to the negative real spectrum.

The scalar objective is `1e-6*sum(sum((spectrum-expt_spec).^2))`, using the transposed experimental data before any display-only processing. For plotting on a non-worker, the experimental matrix is separately denoised with `keep_rank(...,25)`; this cosmetic plot processing does not enter the objective. `fminsearch` is set to iterative display, `MaxIter=10`, `MaxFunEvals=Inf`, and `UseParallel=true`. The source estimates hours while imposing that iteration cap; it does not establish convergence.

## Entry-point output

`glucose_exsy_b()` displays the optimiser's final vector and saves the current figure as `glucose_exsy_b.fig` in the MATLAB working directory. On non-worker evaluations, the two-panel comparison labels the simulated plot “both” and the denoised experimental display “positive”. The function declares no output argument, and the saved figure is not a saved fit-parameter file.

Every objective evaluation rebuilds four directed reaction records from the fitted rates, with corresponding fluorines matched between inside/outside pools. The local rate matrix is retained only for the concentration equilibrium calculation; it is not passed as a retired chemistry input. The initial state uses concentration-weighted `state` without the retired `chem` method.

## Receiver weighting and objective equivalence

The NOESY detection operator must use unweighted `coil_state`; the initial state already carries the equilibrium chemical concentrations. Weighting the receiver again changes the fitted objective. In an integration check with the unweighted production NOESY receiver and the original full acquisition grid, the initial parameter vector gives objective 96212.119538283994, matching the stock concentration convention. This validates the initial objective, not convergence or equivalence of the final optimiser vector; the ten-iteration cap remains unchanged.

This objective requires the unweighted receiver: the reaction-record driver alone does not provide it. A checkout whose `noesy` still constructs detection with `state` weights the concentrations twice and will not reproduce the quoted value; receiver migration is a separate integration prerequisite.
