# examples/fitting/fluoroalkanes/difluoropropane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/difluoropropane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/difluoropropane.m)

- Signature: `difluoropropane()`

## Purpose

Fit 19F and two 1H spectra of 1,3-difluoropropane by varying three scalar-coupling groups and one scale for each spectrum. The source cites [the associated paper](https://doi.org/10.1021/acs.joc.4c00670) and estimates hours of calculation time.

## Experimental data and spin model

The function loads `difluoropropane_fluorine.mat`, `difluoropropane_proton_a.mat`, and `difluoropropane_proton_b.mat`; each contains `axis_ppm` and `spec`. It divides each spectrum by its maximum before fitting. The eight-spin model has six 1H and two 19F spins at field `11.7464`, with S2 groups `[1 2]` and `[7 8]`. Fixed shifts in the source are `4.6075` for protons 1, 2, 7, and 8, `2.1011` for protons 4 and 5, and `-223.5314` for fluorines 3 and 6.

## Fit parameterisation

The initial vector is `[0.3023 1.1605 0.7390 47.0061 5.7870 25.7705]`. Parameters 1–3 scale the 19F, first 1H, and second 1H FIDs; parameters 4–6 set three symmetry-shared groups of scalar couplings. The source uses `fminunc` (unlike the related butane, heptane, and pentane examples, which use `fminsearch`) with `MaxIter=5000` and `MaxFunEvals=Inf`. The objective is the unweighted sum of squared residual norms across the three spectra.

## Acquisition and output

The three non-decoupled liquid acquisitions declare ppm axes and reverse their direction. The 19F sequence sets offset `-105059`, sweep `600`, 2048 points, and zero filling to 4096; 1H A sets offset `2300`, sweep `128`, 512 points, and zero filling to 1024; 1H B sets offset `1050`, sweep `300`, 1024 points, and zero filling to 2048. The source does not label the numeric offset units. The FIDs receive exponential apodisation with source values `8.0`, `9.0`, and `6`; each spectrum is Fourier transformed, reversed to match the axis convention, and interpolated onto its experimental axis with `pchip`. Each objective evaluation plots all three experiment/simulation comparisons and displays the trial vector; the final vector is displayed after optimisation.

The entry point has no output argument, and source code alone supplies no measured best-fit result, uncertainty estimate, or validation claim.
