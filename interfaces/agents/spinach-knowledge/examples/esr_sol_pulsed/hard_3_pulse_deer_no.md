# examples/esr_sol_pulsed/hard_3_pulse_deer_no.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_no.m`
- Signature: `hard_3_pulse_deer_no()`
- Total lines: 56

## Purpose

Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals are positioned at a distance of 25 Angstroms. The numerical calculation is done by brute-force time propaga- tion and numerical powder averaging, including g-factor orien- tation effects on the dipolar coupling. Nitroxide g-tensor is from http://dx.doi.org/10.1063/1.1697233 Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals
- are positioned at a distance of 25 Angstroms.
- The numerical calculation is done by brute-force time propaga-
- tion and numerical powder averaging, including g-factor orien-
- tation effects on the dipolar coupling. Nitroxide g-tensor is
- from http://dx.doi.org/10.1063/1.1697233
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Sequence parameters
- Pulse sequence

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`.
