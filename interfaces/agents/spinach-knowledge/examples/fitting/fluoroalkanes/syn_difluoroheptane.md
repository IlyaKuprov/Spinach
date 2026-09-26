# examples/fitting/fluoroalkanes/syn_difluoroheptane.m

- Signature: `syn_difluoroheptane()`

## Purpose

Fitting of 1H NMR spectrum of syn-3,5-difluoroheptane with respect to J-couplings. See our paper for further details: Methyl groups are ghosted out because they do not influence the signals in question. Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Fitting of 1H NMR spectrum of syn-3,5-difluoroheptane with
- respect to J-couplings. See our paper for further details:
- Methyl groups are ghosted out because they do not influence
- the signals in question.
- Calculation time: hours
- Load experimental data
- Normalize the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result
