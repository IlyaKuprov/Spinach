# kernel/states/stateinfo.m

- Signature: `stateinfo(spin_system,rho,npops)`

## Purpose

Prints the state vector norm and the list of the most populated basis states in the order of decreasing population. Syntax: stateinfo(spin_system,rho,npops)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

## Parameters / inputs

- rho -state vector
- npops -number of largest populations to print

## Outputs

- This function prints a summary of the state composition to the con-
- sole in the following format:
- (L1,M1) (L2,M2) ... (Ln,Mn) coefficient number
- This corresponds to the direct product of single-spin irreducible
- spherical tensors with the specified indices, its coefficient in
- the linear combination, and the number of the corresponding state
- in the basis set.
- Note: this function requires a spherical tensor basis set.

## Implementation structure

- Prints the state vector norm and the list of the most populated basis
- states in the order of decreasing population. Syntax:
- stateinfo(spin_system,rho,npops)
- rho -state vector
- npops -number of largest populations to print
- This function prints a summary of the state composition to the con-
- sole in the following format:
- (L1,M1) (L2,M2) ... (Ln,Mn) coefficient number
- This corresponds to the direct product of single-spin irreducible
- spherical tensors with the specified indices, its coefficient in
- the linear combination, and the number of the corresponding state
- in the basis set.
