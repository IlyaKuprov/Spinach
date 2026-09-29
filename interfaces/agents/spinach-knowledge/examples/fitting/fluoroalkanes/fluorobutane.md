# examples/fitting/fluoroalkanes/fluorobutane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/fluorobutane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/fluorobutane.m)

- Signature: `fluorobutane()`

## Purpose

Fit the source's proton and fluorine spectra by varying scalar couplings and one scale per nucleus. The source header calls the target “2-fluoropentane,” while the function and MAT-file names say “fluorobutane”; the code does not resolve that naming discrepancy. The paper DOI is [10.1021/acs.joc.4c00670](https://doi.org/10.1021/acs.joc.4c00670). The source estimates calculation time in hours.

## Inputs and parameterisation

The entry point loads `fluorobutane_fluorine.mat` (`spec_f`, `axis_f`) and `fluorobutane_proton.mat` (two proton spectra and axes, `spec_ch`/`axis_ch` and `spec_ch2`/`axis_ch2`). Each is normalised by its integral: the first proton segment to −1, the second to −2, and fluorine to −1; the proton segments are concatenated for fitting.

`fminsearch` starts from `[23.9529 6.2219 7.4903 7.1692 4.9608 17.4939 26.2802 48.6869 -14.0924 1.8099 4.3571]`; parameters 1–9 set grouped H–H, H–F, and F–H scalar couplings, and parameters 10 and 11 scale the proton and fluorine simulations. It uses `MaxIter=5000` and unlimited function evaluations. The local objective is the sum of squared residual norms of the real proton and fluorine spectra.

## Spin system and output

The model has nine `1H` spins and one `19F` spin at `sys.magnet=11.7464`, with a Zeeman–Hilbert basis, no approximation, and two `S3` symmetry groups over proton indices 1–3 and 4–6. The proton acquisition uses offset 1400, sweep 2000, 4096 points, and 32768-point zero filling; the separate fluorine acquisition uses offset −81529, sweep 250, 512 points, and 2048-point zero filling. Both use ppm axes, Gaussian apodisation, and `liquid`/`acquire`; simulated axes are interpolated to the experimental axes with `pchip`. Before plotting and computing the objective, the code sets proton-array indices 3770:3820 to zero in both experiment and simulation.

Run `fluorobutane()` with both MAT files available to MATLAB. It has no declared return value: the script displays the final `fminsearch` parameter vector, and the objective displays trial parameters. Its figure compares experiment (red points) and simulation (blue line) in two proton windows (1.5–1.8 and 4.51–4.70 ppm) and a fluorine window (−173.4 to −172.95 ppm).

The source specifies no parameter bounds, uncertainty estimates, or fit-success criterion; the displayed vector is an optimiser result, not a reported measured fit outcome.
