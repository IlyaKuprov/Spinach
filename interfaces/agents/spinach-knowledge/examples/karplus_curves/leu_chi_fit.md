# examples/karplus_curves/leu_chi_fit.m

- Signature: `leu_chi_fit()`
- Source: [examples/karplus_curves/leu_chi_fit.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/karplus_curves/leu_chi_fit.m)

## Method

This driver extracts three Karplus coefficients from a DFT dihedral-angle scan over a leucine chi angle. The source attributes the scan calculations to Gaussian09. It calls karplus_fit('leu_chi_data',{[31 29 23 24]}), passing the named dataset and the atom quartet [31 29 23 24]. The driver does not itself define the fitting function or describe how the data were generated beyond that comment; those details belong to the dataset and karplus_fit implementation.

## Reported quantities

The fitter returns A, B, C and their standard deviations sA, sB, sC. The script prints each coefficient with its standard deviation. No fitted numbers are included here because this driver was not executed and its source contains no hard-coded fit result.

## Scope

The example demonstrates parameter extraction from the supplied scan data, not a new quantum-chemical calculation: this driver only calls the fitter and displays its returned values. It does not state coefficient values, fit quality, uncertainty interpretation, or a DOI.
