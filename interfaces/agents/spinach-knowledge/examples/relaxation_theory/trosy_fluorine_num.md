# examples/relaxation_theory/trosy_fluorine_num.m

- Signature: `trosy_fluorine_num()`

## Purpose

Transverse relaxation rate as a function of the applied magnetic field in a 3-fluorotyrosine labelled protein. The fluorine atom and its directly bonded carbon are included. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Transverse relaxation rate as a function of the applied magnetic
- field in a 3-fluorotyrosine labelled protein. The fluorine atom
- and its directly bonded carbon are included.
- Calculation time: minutes.
- Read 3-fluorotyrosine DFT calculation
- Extract coordinates and CSAs
- Relaxation theory
- Basis set
- Disable startup checks
- Magnetic field grid
- Loop over magnetic fields
- Set the magnet field
