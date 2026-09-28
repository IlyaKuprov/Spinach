# examples/fitting/fluoroalkanes/syn_difluoroheptane.m

- Signature: `syn_difluoroheptane()`

## Purpose

Fit the J-couplings of syn-3,5-difluoroheptane against three experimental spectra: a 19F spectrum (`syn_dfh_fluorine.mat`) and two distinct 1H datasets (`syn_dfh_proton_a.mat` and `syn_dfh_proton_b.mat`). The methyl groups are ghosted out because they do not influence the signals in question. See https://doi.org/doi/10.1021/acs.joc.4c00670 for further details.

Calculation time: hours.

## Workflow

- Load the ppm axes and spectra from all three files and normalize each experimental spectrum by its maximum.
- Start a 15-parameter Nelder–Mead fit with `fminsearch` (`MaxIter=5000`, `MaxFunEvals=Inf`). Parameters 1–3 scale the simulated spectra; parameters 4–15 specify fitted J-couplings, with symmetry-related couplings sharing parameters. Couplings to the ghosted methyl groups are fixed at 7.45.
- For each trial parameter set, build a Spinach spin system at `sys.magnet=11.7464` using the specified isotope, chemical-shift, and coupling assignments; use a Zeeman–Hilbert basis without approximation. Simulate three acquisitions: 19F; 1H A, observing spins 11 and 18; and 1H B, observing spins 8 and 9.
- Apply Gaussian apodisation (15.0, 16.0, and 8.0 for 19F, 1H A, and 1H B, respectively), scale the signals, Fourier-transform and zero-fill them, then interpolate the simulated ppm spectra onto their respective experimental axes.
- Plot experimental and simulated spectra for all three datasets during fitting. Minimize the sum of squared spectral residuals, weighting the 19F residual by 10 and each 1H residual by 1; display the fitted parameter vector.