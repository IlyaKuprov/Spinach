# examples/fitting/fumarate_global.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fumarate_global.m`
- Signature: `fumarate_global()`
- Total lines: 134

## Purpose

Simultaneous fitting of 1H and 13C NMR spectra of a slightly asymmetric fumarate diester. Calculation time: minutes

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Simultaneous fitting of 1H and 13C NMR spectra of a slightly
- asymmetric fumarate diester.
- Calculation time: minutes
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Run the optimisation
- Display the result
- Least squares error function
- Silence Spinach
- Absorb parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `fminsearch()`, `errfun()`, `params()`, `spin()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `sweep2ticks()`, `subplot()`, `kxlabel()`.
