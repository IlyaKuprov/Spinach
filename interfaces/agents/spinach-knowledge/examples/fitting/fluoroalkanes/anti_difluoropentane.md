# examples/fitting/fluoroalkanes/anti_difluoropentane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/anti_difluoropentane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/anti_difluoropentane.m)

- Signature: `anti_difluoropentane()`

## Purpose and naming caveat

The source header describes a fit of the 1H spectrum of **syn**-2,4-difluoropentane, but the callable function and all three experimental-data filenames use **anti**. The source code does not resolve which stereoisomer the measurements represent, so retain both labels rather than silently choosing one. The header cites [the associated paper](https://doi.org/10.1021/acs.joc.4c00670) and estimates hours of calculation time.

## Data and model

The function loads `anti_dfp_fluorine.mat`, `anti_dfp_proton_a.mat`, and `anti_dfp_proton_b.mat`, each containing `axis_ppm` and `spec`. The three spectra are scaled to their maxima with factors `7`, `4`, and `10`, respectively. The model contains ten 1H and two 19F spins at field `11.7464`; its S3 symmetry groups are `[1 2 3]` and `[10 11 12]`. The source fixes the proton shifts at `1.0189` for spins 1–3 and 10–12, `4.9062` for spins 5 and 8, `1.8307` for spins 6 and 7, and the fluorine shifts at `-175.6250` for spins 4 and 9. The basis uses Zeeman-Hilbert formalism without approximation.

## Fit and acquisition workflow

The twelve-entry starting vector is `[6.2155 23.8923 9.8426 2.4324 13.7379 36.4223 49.2708 1.6130 -15.1165 1.0605 0.7229 1.4874]`. Parameters 1–9 set symmetry-related scalar couplings; parameters 10–12 scale the 19F, first 1H, and second 1H FIDs. `fminsearch` minimises the unweighted sum of the three squared spectral-residual norms, with `MaxIter=5000` and `MaxFunEvals=Inf`.

Three non-decoupled acquisitions use the source's ppm-axis and reversed-axis convention. The 19F sequence sets offset `-82646`, sweep `300`, 512 points, and zero filling to 2048; 1H A sets offset `2450`, sweep `180`, 256 points, and zero filling to 1024; 1H B sets offset `915`, sweep `150`, 256 points, and zero filling to 1024. The source does not label the numeric offset units. The FIDs are Gaussian-apodised (source values 7.0, 7.0, and 6.0), Fourier transformed, converted to simulated ppm axes, and interpolated onto the experimental axes with `pchip`. The objective evaluation plots the three experiment/simulation comparisons in stacked panels; the final fit vector is displayed.

The function declares no output argument and the source reports no measured fit outcome, uncertainty, or fit-quality statistic. Numeric offsets are assigned in the acquisition setup, but no offset units are inferred here; the code explicitly sets `axis_units='ppm'`.
