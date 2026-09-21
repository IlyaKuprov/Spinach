# examples/spin_chemistry/cidnp_pumping_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_pumping_1.m`
- Signature: `cidnp_pumping_1()`
- Total lines: 69

## Purpose

A simulation of the matrix in Equation 2 of IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calculation time: seconds.

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A simulation of the matrix in Equation 2 of IK's paper on chemically
- amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
- Calculation time: seconds.
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
- Relaxation theory
- Formalism and basis
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `unit_state()`, `state()`, `magpump()`.
