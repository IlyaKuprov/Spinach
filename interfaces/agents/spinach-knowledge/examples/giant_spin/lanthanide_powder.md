# examples/giant_spin/lanthanide_powder.m

- Signature: `lanthanide_powder()`

## Purpose

Powder spectrum of Gd(III) with ZFS up to 4th spherical rank using the giant spin Hamiltonian formalism in a sweepable 400 MHz NMR magnet and microwaves at 263.2 GHz. Odd ranks are zero because a zero-field Hamiltonian must be even under time reversal. The 4th rank terms are converted from the Stevens parameters b40=4e-4 cm^-1 and b44=-2e-4 cm^-1 reported for Gd(III) in tetragonal BaTiO3 by Rimai and deMars (https://doi.org/10.1103/PhysRev.127.702). Calculation time: seconds.

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Powder spectrum of Gd(III) with ZFS up to 4th spherical rank
- using the giant spin Hamiltonian formalism in a sweepable 400
- MHz NMR magnet and microwaves at 263.2 GHz. Odd ranks are zero
- because a zero-field Hamiltonian must be even under time rever-
- sal. The 4th rank terms are converted from the Stevens parame-
- ters b40=4e-4 cm^-1 and b44=-2e-4 cm^-1 reported for Gd(III)
- in tetragonal BaTiO3 by Rimai and deMars:
- Calculation time: seconds.
- Spin system properties
- Field sweep
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the simulation in the high-T approximation
- Plotting
