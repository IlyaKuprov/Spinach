# examples/fitting/nmr_kinetics/glucose_exsy_a.m

- Signature: `glucose_exsy_a()`

## Purpose

Fits the 2,2,3,3-tetrafluoroglucose NOESY spectrum with respect to chemical-exchange reaction rates and the rotational correlation time in Redfield theory. The source describes a calculation time of minutes and notes that the iteration count is limited in this example.

## Physical and mathematical content

The model combines chemical exchange with Redfield relaxation and rotational motion. A trial parameter vector determines the spin-system and sequence settings used to simulate the NOESY data; the calculated and experimental spectra are compared in a least-squares fit.

## Numerical and algorithmic content

The script optimises an initial guess, simulates the signal, applies apodisation, forms the indirect-dimension States signal, and performs the F2 and F1 Fourier transforms before taking the real spectrum. It processes the experimental spectrum and evaluates the least-squares error. Plotting includes cosmetic SVD denoising and a theory-versus-experiment display.

## Implementation structure

The top-level function sets up the figure, initial guess and optimiser, then calls a local error function. That function configures the Redfield spin system and sequence, runs the simulation, processes both simulated and experimental data, computes and displays the objective and parameters, and plots the comparison.
