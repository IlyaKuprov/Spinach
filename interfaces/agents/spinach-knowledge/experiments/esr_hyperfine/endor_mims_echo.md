# experiments/esr_hyperfine/endor_mims_echo.m

- MATLAB implementation: [experiments/esr_hyperfine/endor_mims_echo.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_hyperfine/endor_mims_echo.m)

Signature: stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)

## Purpose and physical sequence

This routine computes a stimulated-echo diagnostic for the Mims ENDOR sequence with the nuclear RF pulse absent. It prepares electron Lz magnetisation and detects with the electron L+ state. The ideal sequence applies an electron pi/2 (90°) rotation about x, evolves for tau, applies a second x-axis pi/2 rotation, inserts the delay corresponding to the missing nuclear pulse, and applies a negative pi/2 rotation about y. It then records the electron-coil signal over a detection period from 0 to 2*tau. This is a simulated reference echo, not a nuclear-frequency spectrum or measured data.

## Inputs and documented values

- parameters.spins: one-element cell array containing the working spin label, normally {'E'}. For a non-spin-1/2 electron, the source gives {'7E'} as a multiplicity example.
- parameters.electrons: integer indices identifying electron spins in the isotope list; these select the pulse operators, initial state, and detection state.
- parameters.tau: positive delay in seconds between the first two elements of the ENDOR sequence. The earlier page and source document 200e-9 s as the example value.
- parameters.n_dur: positive nuclear-pulse duration, in seconds; 50e-6 s is documented as typical. In this diagnostic it is the duration of the missing nuclear-pulse interval, not an applied RF pulse.
- parameters.nsteps: positive integer number of detection steps over 0 to 2*tau.
- H, R, and K: context-provided Hamiltonian, relaxation, and kinetics matrices with matching dimensions. The function uses Liouville-space formalism and converts to the adjoint representation when needed.

## Output and limitation

stim_echo is the time-sampled stimulated echo in the absence of nuclear RF. The sequence uses ideal rotations and provides a baseline diagnostic; it does not by itself calculate ENDOR modulation, DNP/hyperpolarisation, a magnetic-field sweep, or spatial imaging.

Source: https://spindynamics.org/wiki/index.php?title=endor_mims_echo.m
