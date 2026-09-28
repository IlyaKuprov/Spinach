# examples/karplus_curves/leu_chi_fit.m

- Signature: `leu_chi_fit()`

## Purpose

Fits Karplus coefficients to the DFT dihedral-angle scan for a leucine chi angle, using Gaussian09-derived data. The script calls `karplus_fit('leu_chi_data',{[31 29 23 24]})` and displays fitted `A`, `B`, and `C` coefficients with their standard deviations.

## Physical / mathematical content

A Karplus fit relates torsion angle to scalar coupling; this example extracts the three coefficients and their uncertainties from the supplied scan data.

## Implementation structure

Calls `karplus_fit` for the `leu_chi_data` dataset and atom quartet `[31 29 23 24]`, then prints each coefficient and standard deviation.
