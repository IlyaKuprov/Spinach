# examples/fitting/maleate_global.m

- Signature: `maleate_global()`

## Purpose

Simultaneously fits the 1H and 13C NMR spectra of a slightly asymmetric maleate diester. The source estimates a calculation time of hours.

## Physical and mathematical content

A shared parameter vector defines the spin-system model used to fit both experimental nuclei. The source loads and normalises the proton and carbon data separately, then evaluates both simulated spectra in a joint least-squares objective.

## Numerical and algorithmic content

The top-level function supplies an initial guess and optimiser options, minimises the combined spectral residual, and displays the fitted parameters. In the local error function, each trial vector is absorbed into a Spinach system; the two acquisitions are simulated, Fourier processed and aligned with their measured axes before their residuals are combined.

## Implementation structure

The workflow proceeds from loading the two spectra to optimisation. The error function configures the spin system, basis and separate 1H/13C acquisition parameters, computes and processes both signals, plots theory against experiment, and returns the objective value.
