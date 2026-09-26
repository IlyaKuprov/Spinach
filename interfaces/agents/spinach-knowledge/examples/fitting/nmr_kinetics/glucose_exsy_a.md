# examples/fitting/nmr_kinetics/glucose_exsy_a.m

- Signature: `glucose_exsy_a()`

## Purpose

Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the reaction rates in a chemical exchange problem and the rotational correlation time within Redfield theory. Calculation time: minutes (iteration count is limited in this example file)

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the
- reaction rates in a chemical exchange problem and the rotational
- correlation time within Redfield theory.
- Calculation time: minutes (iteration count is limited
- in this example file)
- Get a figure going
- Set the initial guess
- Set optimiser options
- Run the optimisation
- Display the result
- Hush up Spinach
- Magnet field
