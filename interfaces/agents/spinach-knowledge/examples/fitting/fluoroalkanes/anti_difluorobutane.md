# examples/fitting/fluoroalkanes/anti_difluorobutane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/anti_difluorobutane.m`
- Signature: `anti_difluorobutane()`
- Total lines: 158

## Purpose

Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise the data
- Concatenate spectral intervals
- Broaden out Z1 shim problem
- Set the guess
- Get a figure going
- Set optimiser options
- Run the optimisation
- Display the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `trapz()`, `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `sim_spec()`, `sweep2ticks()`.
