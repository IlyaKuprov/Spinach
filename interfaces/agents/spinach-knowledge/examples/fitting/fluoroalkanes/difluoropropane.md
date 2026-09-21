# examples/fitting/fluoroalkanes/difluoropropane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/difluoropropane.m`
- Signature: `difluoropropane()`
- Total lines: 175

## Purpose

Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result
- Least squares error function
- Silence the output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminunc()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spec_theo_f()`, `spec_theo_ha()`, `spec_theo_hb()`.
