# examples/fitting/fluoroalkanes/syn_difluoropentane.m

- Signature: `syn_difluoropentane()`

## Purpose

Fits the 1H NMR spectrum of syn-2,4-difluoropentane with respect to J-couplings. The source also simulates and fits the corresponding 19F and two 1H data sets. See the paper: https://doi.org/10.1021/acs.joc.4c00670. The source notes a calculation time of hours.

## Physical and mathematical content

The model contains ten 1H and two 19F spins, with chemical shifts and scalar couplings parameterised for the fit. The two equivalent three-proton groups are represented with S3 symmetry. Experimental spectra are loaded, scaled/shifted, and compared with Spinach simulations; the fitted vector controls couplings and signal scaling.

## Numerical and algorithmic content

The objective is the sum of squared spectral residual norms for the 19F spectrum and both 1H spectra. The script searches the supplied initial parameter vector with `fminsearch`. Simulated FIDs are apodised, Fourier transformed, converted to frequency axes, interpolated onto the experimental axes, and plotted against the data.

## Implementation structure

The function loads the three experimental data files, prepares the initial guess and optimiser options, then calls the local error function. That function builds the spin system and symmetry-adapted basis, configures separate 19F and 1H acquisitions, simulates the three spectra, processes them, and returns the combined least-squares error.
