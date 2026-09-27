# examples/fitting/fluoroalkanes/anti_difluoropentane.m

- Signature: `anti_difluoropentane()`

## Purpose

Fit simulated ¹⁹F and two ¹H NMR spectral regions to experimental data by varying J-couplings and three spectral amplitudes. The source header describes this as fitting the ¹H spectrum of **syn**-2,4-difluoropentane, whereas the function and data files are named **anti**; the code does not establish which stereoisomer the data represent. See the paper: https://doi.org/doi/10.1021/acs.joc.4c00670

Calculation time: hours.

## Workflow

1. Load `axis_ppm` and `spec` from `anti_dfp_fluorine.mat`, `anti_dfp_proton_a.mat`, and `anti_dfp_proton_b.mat`. Scale the three experimental spectra by `7/max(spec)`, `4/max(spec)`, and `10/max(spec)`, respectively.
2. Start Nelder–Mead optimisation with `fminsearch` from `[6.2155, 23.8923, 9.8426, 2.4324, 13.7379, 36.4223, 49.2708, 1.6130, -15.1165, 1.0605, 0.7229, 1.4874]`. Optimiser settings are `Display='iter'`, `MaxIter=5000`, and `MaxFunEvals=Inf`; the fitted vector is displayed.
3. For each trial vector, construct a 12-spin system containing ten ¹H and two ¹⁹F spins at a magnetic field of `11.7464`. Fixed chemical shifts are `1.0189` for spins 1–3 and 10–12, `4.9062` for spins 5 and 8, `1.8307` for spins 6 and 7, and `-175.6250` for spins 4 and 9. Parameters 1–9 set symmetry-related J-couplings; the coupling between spins 4 and 9 is parameter 8, and that between spins 6 and 7 is parameter 9. Use the `zeeman-hilb` formalism without approximation, with `S3` symmetry for spin groups `[1 2 3]` and `[10 11 12]`.
4. Simulate three acquisitions without decoupling: ¹⁹F excitation/detection over the ¹⁹F spins (`offset=-82646`, `sweep=300`, `npoints=512`, `zerofill=2048`); ¹H A excitation/detection on spins `[5 8]` (`offset=2450`, `sweep=180`, `npoints=256`, `zerofill=1024`); and ¹H B on `[6 7]` (`offset=915`, `sweep=150`, `npoints=256`, `zerofill=1024`). All three specify ppm axes and axis inversion.
5. Apply Gaussian apodisation of `7.0`, `7.0`, and `6.0`, respectively, multiplying by fitted amplitude parameters 10–12 and dividing by `4e3`. Take the real, shifted, zero-filled FFT and reverse each spectrum. Convert simulated frequency ticks to ppm using the corresponding spin base frequency, then interpolate onto each experimental axis with `pchip`.
6. Plot experimental points against simulated curves in three stacked, reverse-ppm-axis panels. Minimise the sum of squared Euclidean residual norms across all three spectra.