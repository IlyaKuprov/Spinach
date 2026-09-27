# examples/fitting/nmr_kinetics/glucose_exsy_b.m

- Signature: `glucose_exsy_b()`

## Purpose

Fits the 3,3-difluoroglucose NOESY spectrum with respect to chemical-exchange reaction rates and rotational correlation times in Redfield theory. The source estimates a calculation time of hours and notes that the iteration count is limited in this example.

## Physical and mathematical content

The fit treats chemical exchange and Redfield relaxation, with rotational correlation times among the varied parameters. A trial parameter vector sets up the spin system and NOESY simulation; simulated and experimental spectra are compared using a least-squares objective.

## Numerical and algorithmic content

The simulated signal is apodised and processed by F2 and F1 Fourier transforms using the States signal, after which the real spectrum is compared with the processed experimental spectrum. The script reports the error and parameters and plots theory against experiment; its figure path includes cosmetic SVD denoising.

## Implementation structure

The function sets the figure, initial guess and optimiser options, then runs the fit through a local error function. The error function configures the Redfield spin-system and sequence parameters, performs the simulation and spectral processing, evaluates the least-squares mismatch, and produces the comparison plot.
