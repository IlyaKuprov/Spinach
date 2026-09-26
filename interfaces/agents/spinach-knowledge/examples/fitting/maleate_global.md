# examples/fitting/maleate_global.m

- Signature: `maleate_global()`

## Purpose

Simultaneous fitting of 1H and 13C NMR spectra of a slightly asymmetric maleate diester. Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Simultaneous fitting of 1H and 13C NMR spectra of a slightly
- asymmetric maleate diester.
- Calculation time: hours
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Run the optimisation
- Display the result
- Least squares error function
- Silence Spinach
- Absorb parameters
