# examples/fitting/fluoroalkanes/syn_difluoropentane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/syn_difluoropentane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/syn_difluoropentane.m)

- Signature: `syn_difluoropentane()`

## Purpose and experimental convention

Fit the three source-provided spectra for syn-2,4-difluoropentane: a 19F trace and two 1H traces. The source scales the experimental traces before fitting: fluorine to peak height 7, proton A to peak height 4 followed by a −0.1 offset, and proton B to peak height 10. See [10.1021/acs.joc.4c00670](https://doi.org/10.1021/acs.joc.4c00670); the source estimates hours of calculation.

The entry point loads `syn_dfp_fluorine.mat`, `syn_dfp_proton_a.mat`, and `syn_dfp_proton_b.mat`, each with `axis_ppm` and `spec`. The source initial guess is `[6.2780 23.5067 5.1361 7.0537 24.9899 16.9264 48.0475 1.6017 -14.5413 0.7437 0.8771 0.6311]`; `fminsearch` uses `MaxIter=5000` and unlimited function evaluations. The first nine parameters populate grouped scalar couplings; the final three scale the 19F, 1H A, and 1H B simulations.

## Spin system and simulation

The model has ten `1H` and two `19F` spins at `sys.magnet=11.7464`, with a Zeeman–Hilbert basis, no approximation, and `S3` symmetry for proton indices 1–3 and 10–12. Three `liquid`/`acquire` simulations observe 19F, 1H A (spins 4 and 8), and 1H B (spins 6 and 7), with no decoupling. The source sets their offset/sweep/point-count/zero-fill tuples to (−81655, 300, 512, 2048), (2426, 128, 256, 1024), and (986, 350, 512, 2048), respectively; axes are labelled ppm. Fixed Gaussian apodisation arguments are 7, 7, and 6 for the channels, with fitted scale factors applied before Fourier transformation. The theoretical spectra are converted to ppm and interpolated to their experimental axes with `pchip`.

## Entry point and output

Run `syn_difluoropentane()` with all three MAT files available. The objective is the unweighted sum of squared residual norms over the three spectra. The function returns no declared output; it displays the parameter vector and produces three reversed-ppm panels comparing experiment and simulation. No fitted numerical outcome, uncertainty estimate, or success criterion is given in the source. The source does not attach units to its Gaussian apodisation arguments.
