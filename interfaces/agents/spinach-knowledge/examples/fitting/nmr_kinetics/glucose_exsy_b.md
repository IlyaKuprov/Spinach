# examples/fitting/nmr_kinetics/glucose_exsy_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/nmr_kinetics/glucose_exsy_b.m`
- Signature: `glucose_exsy_b()`
- Total lines: 172

## Purpose

Fitting of 3,3-difluoroglucose NOESY with respect to the reaction rates in a chemical exchange and the rotational correlation times within Redfield theory. Calculation time: hours (iteration count is limited in this example file)

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `exsy_err()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Fitting of 3,3-difluoroglucose NOESY with respect to the reaction
- rates in a chemical exchange and the rotational correlation times
- within Redfield theory.
- Calculation time: hours (iteration count is limited
- in this example file)
- Get a figure going
- Set the initial guess
- Set optimiser options
- Run the optimisation
- Display the result
- Save figure
- Hush up Spinach

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `savefig()`, `exsy_err()`, `num2cell()`, `params()`, `equilibrate()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`.
