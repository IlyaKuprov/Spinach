# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m

- Signature: `mas_powder_mhc_fplanck_exchange()`
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m

## Purpose

Simulates the MAS 1H NMR spectrum of water protons in monohydrocalcite with position exchange between two reaction endpoints. The source cites https://doi.org/10.1038/s41467-023-44381-x, attributes its parametrisation to Huang et al. (ACIE, 2021), and reports seconds of calculation time.

## Spin and exchange model

The code reads the CASTEP-derived `mhc.magres` file and removes C, O, and Ca. It builds two two-proton endpoints from proton sites at positions 1 and 4, with the site order reversed in the second endpoint; all four spins are 1H. Their shielding tensors are `29.25*eye(3)-cst` using the corresponding source tensors, and the coordinates are selected from those two positions. The kinetic parts are `[1 2]` and `[3 4]`, with initial concentrations `[1 1]` in arbitrary units. The exchange-rate matrix is `2,000*[-1 1; 1 -1]` Hz. The field setting is `sys.magnet=9.4`, labelled 400 MHz NMR in the source; the basis is `sphten-liouv` with no approximation.

## MAS and pulse-acquire settings

The source sets the MAS rate parameter to 10,000 (no unit is written beside this assignment) about axis `[1 1 1]`, with powder grid `rep_2ang_100pts_sph` and maximum rank 13. The sweep is set as `1/(5e-6)`; the source does not annotate a unit on this assignment. It uses 512 points and 1,024-point zero filling. The initial state and detection operator are both `state(spin_system,'L+','1H')`, and the sequence callback is `@acquire` inside `singlerot`. No RF pulse duration or power is specified in this source.

## Inputs and outputs

The simulation applies exponential apodisation with parameter 6 to the computed FID, Fourier transforms it, and plots the real spectrum with `plot_1d`. Its input is the CASTEP-derived structural/shielding data in `mhc.magres`; it does not load an experimental FID or measured spectrum. The output is a simulated MAS spectrum.
