# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m

- Signature: `mas_powder_mhc_fplanck()`
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m

## Purpose

Simulates the magic-angle-spinning 1H NMR spectrum of all protons in the monohydrocalcite unit cell. The source cites https://doi.org/10.1038/s41467-023-44381-x and estimates hours of runtime, or minutes on a GPU.

## Spin model

The calculation reads the CASTEP-derived `mhc.magres` file, removes C, O, and Ca, and assigns the remaining 18 sites as 1H. It converts each proton shielding tensor with the Huang et al. ACIE 2021 parametrisation `29.25*eye(3)-cst` and supplies the selected Cartesian coordinates. The source sets `sys.magnet=9.4` and labels this as 400 MHz NMR. The basis is `sphten-liouv` with approximation `IK-0`, interaction level 3, and projection `+1`; the interaction cutoff is 500 Hz.

## MAS and pulse-acquire settings

The source sets the MAS rate parameter to 10,000 (no unit is written beside this assignment) and the rotor axis to `[1 1 1]`. The powder grid is `rep_2ang_100pts_sph`, with maximum rank 13. The acquisition sweep is set as `1/(5e-6)`; the source does not annotate a unit on this assignment. It uses 512 points and 1,024-point zero filling. Both the initial state and detection operator are `state(spin_system,'L+','1H')`. The source uses the pulse-acquire callback `@acquire` through `singlerot`; it does not specify an RF pulse duration or power in this file.

## Inputs and outputs

The code simulates an FID, applies exponential apodisation with parameter 6, Fourier transforms it, and plots the real spectrum with `plot_1d`. The input is the CASTEP-derived `mhc.magres` data; no experimental FID or measured spectrum is loaded here. The output is the simulated MAS spectrum. The source's GPU-enable line is commented out.
