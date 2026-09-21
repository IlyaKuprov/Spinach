# kernel/states/stateinfo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/stateinfo.m`
- Signature: `stateinfo(spin_system,rho,npops)`
- Total lines: 83

## Purpose

Prints the state vector norm and the list of the most populated basis states in the order of decreasing population. Syntax: stateinfo(spin_system,rho,npops)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `rho()`, `sorting_index()`, `states()`, `lin2lm()`, `largest_states()`, `cell2mat()`, `largest_elemts()`, `ismember()`, `iscolumn()`, `isscalar()`.
