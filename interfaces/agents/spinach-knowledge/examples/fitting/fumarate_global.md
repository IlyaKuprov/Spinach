# examples/fitting/fumarate_global.m

- Signature: `fumarate_global()`

## Purpose

Simultaneously fits the 1H and 13C NMR spectra of a slightly asymmetric fumarate diester. The source estimates a calculation time of minutes.

## Physical and mathematical content

The fitted model uses one parameter vector to describe the fumarate spin system for both nuclei. The experimental proton and carbon spectra are loaded and normalised; Spinach simulations for each channel are compared with their corresponding data through a joint least-squares objective.

## Numerical and algorithmic content

The script sets an initial guess and optimiser options, minimises the spectral residual objective, and displays the resulting parameter vector. The error function constructs the spin system, simulates the 1H and 13C acquisitions, processes the FIDs into spectra, and combines the two residual contributions.

## Implementation structure

The top-level function loads the two data sets and runs the optimisation. Its local error function absorbs the trial parameters, defines the spin-system and acquisition settings for each nucleus, performs both simulations, transforms and aligns their frequency axes with the experimental axes, plots experiment against simulation, and returns the least-squares error.
