# examples/karplus_curves/tyr_chi_fit.m

## Purpose and callable context

The no-argument MATLAB entry point `tyr_chi_fit()` fits Karplus coefficients for a DFT dihedral-angle scan over one tyrosine chi angle. The source comments identify Gaussian09 as the calculation package used for the scan. It calls `karplus_fit('.\tyr_chi_data',{[15 14 11 12]})`: the relative data directory is `tyr_chi_data` and the four atom indices define the fitted torsion. Running the wrapper therefore depends on the `karplus_fit` helper and that dataset being available in the execution context.

## Method and output

The helper returns `A`, `B`, `C` and their reported standard deviations `sA`, `sB`, `sC`; the wrapper prints those six values. The source does not state the fit equation or the units, and contains no fitted numerical results. It makes no plot.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/karplus_curves/tyr_chi_fit.m)
