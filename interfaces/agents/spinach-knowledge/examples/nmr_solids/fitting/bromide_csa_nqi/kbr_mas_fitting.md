# examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m

- Signature: kbr_mas_fitting()

## Purpose and scope

This example fits a ⁷⁹Br MAS NMR spectrum of potassium bromide to quadrupolar-coupling parameters. The source says a single tensor does not fit the spectrum and that at least three are necessary, likely reflecting a distribution of electrostatic environments in the powder. It estimates hours of calculation. These are source comments, not a fit or validation performed here.

## Spin model and code-set inputs

The input file is KBr_400MHz_2kHz.txt; column 2 is flipped and divided by 100 to form the comparison trace. The source sets the magnet field to 9.3659 T, MAS rate to 2000 Hz, receiver offset to 6034.96 Hz, sweep width to 1e5 Hz, 4096 points, and zero-fill to 32768.

The model contains three ⁷⁹Br spins with a shared isotropic chemical shift and three diagonal quadrupolar tensors. For each tensor the third diagonal element is minus the sum of the first two, and the initial guess supplies the XX and YY values in kHz before the code multiplies them by 1e3. The basis is sphten-liouv with IK-0 approximation, inter-level 1 and projection +1. One fitted relaxation-rate value is assigned to both T1 and T2 for all three spins; the model uses zero equilibrium and diagonal relaxation retention. The starting density operator is a weighted sum of the three L+ states, with the three component weights among the fit parameters.

The initial parameter vector is [60.0933, 13.7569, 1.6424, 4.0779, 4.5179, 1.5885, 0.9449, 263.9835, 40, 32, 28]: the source labels the first value as chemical shift in ppm, six tensor eigenvalue inputs in kHz, the relaxation rate in Hz, and the last three as weights. They are optimiser starting values, not reported fitted values.

## MAS calculation, comparison, and limits

For each trial parameter vector the script simulates a single-rotation spectrum using singlerot and acquire, with the rep_2ang_200pts_oct powder grid, maximum rank 50, 79Br observation, and axis units set to Hz. It plots the normalised input trace as red points and the calculated real Fourier spectrum as a blue line. The least-squares objective is the squared 2-norm of their difference; fminsearch is configured for up to 5000 iterations.

The code models a one-dimensional MAS spectrum; it does not define a CP, HMQC, or Hartmann–Hahn transfer block. The wrapper call supplies the acquisition callback and spectral settings, not an experimental pulse sequence beyond what is explicit in the source. The file contains no saved best-fit parameter set or measured output. The source file gives no DOI.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/fitting/bromide_csa_nqi/kbr_mas_fitting.m
