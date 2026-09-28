# examples/optimal_control/distortions/distortions_figure_3.m

- Signature: `distortions_figure_3()`

## Purpose

Figure 3 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- The script optimises a 45-slice, two-channel `13C` GRAPE pulse using `fmaxnewton` with `@grape_xy` and the `lbfgs` method. The last five slices are frozen as dead time.
- It compares an undistorted GRAPE pulse, a conventional square pulse, a GRAPE pulse distorted by two successive RLC filters, and a GRAPE pulse optimised with those filters included in the control model. The plotted comparison is labelled “RAW-GRAPE (RLC distorted)”.

## Numerical / algorithmic content

- Each simulated free induction decay is Gaussian-apodised, Fourier-transformed with zero filling to 16,384 points, and plotted as the real spectrum; the acquisition sweep is 70,000 with 2,048 points.

## Implementation structure

- Figure 3 from the paper by Rasulov and Kuprov:
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping
- Set up spin states
- Get the control operators
- Get the drift Hamiltonian
- Define control parameters
- Last five slices are dead time
