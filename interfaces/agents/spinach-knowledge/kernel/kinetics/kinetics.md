# kernel/kinetics/kinetics.m

- Signature: `K=kinetics(spin_system)`

## Purpose

Chemical kinetics superoperator. Syntax: K=kinetics(spin_system)

## Physical / mathematical content

- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- spin_system -Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual. All adjustable pa-
- rameters are described in the chemical
- kinetics parameters section.

## Outputs

- K -kinetics superoperator. If a Liouvillian is
- assembled manually, this dissipative super-
- operator must enter as 1i*K, for example
- L=H+1i*R+1i*K
- Note: a large variety of chemical reaction models is supported,
- see the chemical kinetics parameters section of the onli-
- ne manual.
- Note: Spinach context functions include relaxation and kinetics
- superoperators into the total Liovillian automatically.

## Implementation structure

- Chemical kinetics superoperator. Syntax:
- K=kinetics(spin_system)
- spin_system - Spinach spin system description object
- produced as described in the spin system
- and basis specification sections, of the
- of the online manual. All adjustable pa-
- rameters are described in the chemical
- kinetics parameters section.
- K - kinetics superoperator. If a Liouvillian is
- assembled manually, this dissipative super-
- operator must enter as 1i*K, for example
- L=H+1i*R+1i*K
