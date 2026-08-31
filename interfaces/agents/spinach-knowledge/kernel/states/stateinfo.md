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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(spin_system,rho,npops)`.
- Lines 36-37: Report the state vector norm; implemented by `report(spin_system,['state vector 2-norm: ' num2str(norm(rho,2))])`.
- Lines 39-40: Locate npops most populated states and sort by amplitude; implemented by `[~,sorting_index]=sort(abs(rho),1,'descend')`.
- Lines 44-45: Print the states and their populations; implemented by `report(spin_system,[num2str(npops) ' most populated basis states (state, coeff, number)'])`.

### Control flow inferred from the code

- Line 46: `for` loop over `n=1:npops`.
- Line 48: `for` loop over `k=1:spin_system.comp.nspins`.
- Line 50: conditional branch on `l==0`.

### Key state/data transformations

- Lines 40: computes `[~,sorting_index]` using `[~,sorting_index]=sort(abs(rho),1,'descend')`.
- Lines 41: computes `largest_elemts` using `largest_elemts=rho(sorting_index(1:npops))`.
- Lines 42: computes `largest_states` using `largest_states=spin_system.bas.basis(sorting_index(1:npops),:)`.
- Lines 47: computes `state_string` using `state_string=cell(1,spin_system.comp.nspins)`.
- Lines 49: computes `[l,m]` using `[l,m]=lin2lm(largest_states(n,k))`.
- Lines 51: computes `state_string{k}` using `state_string{k}=' . '`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(spin_system,rho,npops)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('state analysis is only available for sphten-liouv formalism.')`.

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
