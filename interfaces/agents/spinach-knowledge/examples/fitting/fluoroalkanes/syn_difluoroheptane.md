# examples/fitting/fluoroalkanes/syn_difluoroheptane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/syn_difluoroheptane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/syn_difluoroheptane.m)

- Signature: `syn_difluoroheptane()`

## Purpose and data

Fit the spectra assigned by the source to syn-3,5-difluoroheptane: one 19F trace and two distinct 1H traces. The source comment says methyl groups are ghosted because they do not influence the signals in question. See [10.1021/acs.joc.4c00670](https://doi.org/10.1021/acs.joc.4c00670); its runtime estimate is hours.

The entry point loads `syn_dfh_fluorine.mat`, `syn_dfh_proton_a.mat`, and `syn_dfh_proton_b.mat` (each supplies `axis_ppm` and `spec`) and normalises each experimental spectrum by its maximum. The 15-element initial guess `[1.0552 1.1935 0.9043 -13.8188 4.8881 7.0776 4.4430 7.7255 -14.7896 18.3289 25.8358 48.5693 17.2030 30.4926 1.8459]` feeds `fminsearch` (`MaxIter=5000`, unlimited function evaluations): parameters 1–3 scale the three simulated data sets, while parameters 4–15 provide fitted scalar couplings, with equal parameters reused across paired couplings. Couplings to ghost spins are assigned 7.45 in the source rather than fitted.

## Spin system and acquisitions

The model uses `sys.magnet=11.7464`, isotope labels including `G` ghost spins, the source-assigned chemical shifts, and a Zeeman–Hilbert basis with no approximation. Separate liquid-state acquisitions simulate 19F, 1H A (observed spins 11 and 18), and 1H B (observed spins 8 and 9). Their source parameters are, respectively: offsets −85656, 2300, and 950; sweeps 300, 128, and 300; point counts 512, 256, and 512; and zero fills 2048, 1024, and 2048. Axes are labelled ppm. Each FID receives a fixed Gaussian apodisation argument (15, 16, or 8) as well as its fitted scale, then is Fourier transformed and interpolated onto its experimental axis with `pchip`.

## Entry point and objective

Run `syn_difluoroheptane()` with the three MAT files available. The function displays the current and final parameter vectors and plots experiment against simulation in three reversed-ppm panels. Its objective is the sum of squared residual norms, weighting 19F by 10 and each 1H trace by 1; it has no declared return value. The source reports no fit result, uncertainty, or success threshold, and does not state units for the Gaussian apodisation arguments.
