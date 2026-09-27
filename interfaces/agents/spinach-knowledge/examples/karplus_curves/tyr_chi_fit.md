# examples/karplus_curves/tyr_chi_fit.m

- Signature: `tyr_chi_fit()`

## Purpose

Fits Karplus coefficients to the DFT dihedral-angle scan for a tyrosine chi angle, using Gaussian09-derived data. The script calls `karplus_fit('.\tyr_chi_data',{[15 14 11 12]})` and displays fitted `A`, `B`, and `C` coefficients with their standard deviations.

## Physical / mathematical content

A Karplus fit relates torsion angle to scalar coupling; this example extracts the three coefficients and their uncertainties from the supplied scan data.

## Implementation structure

Calls `karplus_fit` for the `tyr_chi_data` dataset and atom quartet `[15 14 11 12]`, then prints each coefficient and standard deviation.
