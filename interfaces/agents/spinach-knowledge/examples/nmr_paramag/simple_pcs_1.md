# examples/nmr_paramag/simple_pcs_1.m

- Signature: `simple_pcs_1()`

## Purpose

Pseudocontact shift and Curie relaxation on a proton due to the presence of a point magnetic susceptibility centre. Calculation time: seconds

## Physical / mathematical content

- Paramagnetic NMR examples. These files work with pseudocontact shifts, paramagnetic relaxation, susceptibility tensors, and inverse problems for metal-site localisation or distributed electron density reconstruction.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Pseudocontact shift and Curie relaxation on a proton due to
- the presence of a point magnetic susceptibility centre.
- Calculation time: seconds
- System specification
- Diamagnetic shifts and coordinates
- Magnetic susceptibility
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Spinach relaxation rates
- Summary
