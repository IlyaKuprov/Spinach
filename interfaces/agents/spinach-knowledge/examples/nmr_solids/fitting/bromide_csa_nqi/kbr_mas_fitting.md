# examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m

- Signature: `kbr_mas_fitting()`

## Purpose

Fitting of a 79Br MAS NMR spectrum of potassium bromide with respect to the quadrupole coupling constant. The spectrum cannot be fitted with a single quadrupolar tensor; at least 3 are necessary, likely due to a distribution of electrostatic environments in the powder. Calculation time: hours.

## Physical / mathematical content

The fitting model contains three 79Br sites with a shared isotropic chemical shift and three diagonal quadrupolar coupling tensors. The source comments that one tensor cannot fit the spectrum and that at least three are needed, likely because the powder has a distribution of electrostatic environments.

## Numerical / algorithmic content

The script reads `KBr_400MHz_2kHz.txt`, uses a 2 kHz MAS rate, and sets the acquisition sweep width to 100 kHz (`parameters.sweep=1e5 Hz`). It optimises the model parameters with `fminsearch` (up to 5000 iterations). It simulates with an IK-0 basis, rank 50 and `rep_2ang_200pts_oct`; the objective is the squared 2-norm between experimental and calculated spectra. The source disables hygiene and trajectory-level output inside the objective and plots the experimental and fitted curves at each evaluation.

## Implementation structure

Loads and normalises the experimental spectrum, sets the instrumental parameters and initial guesses, then repeatedly constructs the three-site spin system, simulates and Fourier transforms its signal, and evaluates the least-squares residual.
