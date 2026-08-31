# kernel/reduce.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/reduce.m`
- Signature: `projectors=reduce(spin_system,L,rho)`
- Total lines: 315

## Purpose

Symmetry and trajectory-level state space reduction. Tries all applicable reduction methods (unless disabled during the call to create.m) and returns a cell array of projectors into a set of independently evolving reduced subspaces. Syntax: projectors=reduce(spin_system,L,rho)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check the input; implemented by `grumble(spin_system,L,rho)`.
- Lines 46-47: Check for blanket ban; implemented by `if ismember('trajlevel',spin_system.sys.disable)`.
- Lines 52-53: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 57-58: If a cell array is supplied, make a representative density matrix; implemented by `if iscell(rho)`.
- Lines 66-67: Run symmetry factorization; implemented by `if ismember('symmetry',spin_system.sys.disable)`.
- Lines 69-70: Issue a reminder to the user; implemented by `report(spin_system,'WARNING - permutation symmetry treatment disabled by user.')`.
- Lines 72-73: Return unit matrix; implemented by `projectors{1}=1`.
- Lines 77-78: Inform the user; implemented by `report(spin_system,'no permutation symmetry information has been supplied.')`.
- Lines 85-86: Decide which irreps to keep; implemented by `n_irreps=numel(spin_system.bas.irrep)`.
- Lines 90-91: Check the irrep contribution to the total norm; implemented by `if spin_system.bas.irrep(n).dimension==0`.
- Lines 93-95: Update the user; implemented by `report(spin_system,['irrep #' num2str(n) ' has no states in it - dropped.'])`.
- Lines 97-98: Flag the irrep for dropping; implemented by `irrep_keep_index(n)=0`.
- Lines 103-106: Update the user; implemented by `report(spin_system,['irrep #' num2str(n) ' has less than ' num2str(spin_system.tols.irrep_drop) ' of the initial state norm - dropped.'])`.
- Lines 113-116: Update the user; implemented by `report(spin_system,['irrep #' num2str(n) ' contains ' num2str(spin_system.bas.irrep(n).dimension) ' states - kept.'])`.
- Lines 122-123: Compile the projector array; implemented by `projectors={spin_system.bas.irrep(irrep_keep_index).projector}`.
- Lines 129-130: If a stack is supplied, choose a representative wavefunction; implemented by `if size(rho,2)>1, rho=mean(abs(rho),2); end`.
- Lines 159-160: Update the user; implemented by `report(spin_system,['irrep #' num2str(n) ' has dimension zero - dropped.'])`.
- Lines 167-170: Update the user; implemented by `report(spin_system,['irrep #' num2str(n) ', dimension ' num2str(spin_system.bas.irrep(n).dimension) ' has less than ' num2str(spin_system.tols.irrep_drop) ' of the state…`.

### Control flow inferred from the code

- Line 47: conditional branch on `ismember('trajlevel',spin_system.sys.disable)`.
- Line 53: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `'zeeman-wavef'`, `{'zeeman-liouv','sphten-liouv'}`.
- Line 58: conditional branch on `iscell(rho)`.
- Line 60: `for` loop over `n=2:numel(rho)`.
- Line 67: conditional branch on `ismember('symmetry',spin_system.sys.disable)`.
- Line 88: `for` loop over `n=1:n_irreps`.
- Line 91: conditional branch on `spin_system.bas.irrep(n).dimension==0`.
- Line 130: conditional branch on `size(rho,2)>1, rho=mean(abs(rho),2); end`.
- Line 133: conditional branch on `ismember('symmetry',spin_system.sys.disable)`.
- Line 154: `for` loop over `n=1:n_irreps`.
- Line 157: conditional branch on `spin_system.bas.irrep(n).dimension==0`.
- Line 193: conditional branch on `size(rho,2)>1, rho=mean(abs(rho),2); end`.
- Line 196: conditional branch on `ismember('symmetry',spin_system.sys.disable)`.
- Line 218: `for` loop over `n=1:n_irreps`.

### Key state/data transformations

- Lines 49: computes `projectors{1}` using `projectors{1}=1; return`.
- Lines 59: computes `rho_rep` using `rho_rep=abs(rho{1})`.
- Lines 63: computes `rho` using `rho=rho_rep/numel(rho)`.
- Lines 86: computes `n_irreps` using `n_irreps=numel(spin_system.bas.irrep)`.
- Lines 87: computes `irrep_keep_index` using `irrep_keep_index=true(n_irreps,1)`.
- Lines 98: computes `irrep_keep_index(n)` using `irrep_keep_index(n)=0`.
- Lines 123: computes `projectors` using `projectors={spin_system.bas.irrep(irrep_keep_index).projector}`.
- Lines 259: computes `zte_projector` using `zte_projector=zte(spin_system,projectors{n}'*L*projectors{n},projectors{n}'*rho)`.
- Lines 262: computes `projectors{n}` using `projectors{n}=projectors{n}*zte_projector`.
- Lines 275: computes `pt_projectors` using `pt_projectors=path_trace(spin_system,projectors{n}'*L*projectors{n},projectors{n}'*rho)`.
- Lines 279: computes `pt_projectors{k}` using `pt_projectors{k}=projectors{n}*pt_projectors{k}`.

### Local helper functions

- Line 298: `grumble()` — `function grumble(spin_system,L,rho)`. For centuries, the battle of morality was fought between those who claimed that your life belongs to God and those who claimed that it belongs to your
  - Representative operation: `if ~isnumeric(L)`.
  - Representative operation: `error('L must be numeric.')`.

## Parameters / inputs

- L -Liouvillian matrix
- rho -initial state (source state screening) or
- destination state (destination state screening)

## Outputs

- projectors -a cell array of projectors into independently
- evolving reduced subspaces. The projectors are
- to be used as follows:
- L_reduced=P'*L*P; (for matrices)
- rho_reduced=P'*rho; (for state vectors)
- Notes: further information on what this function does is avai-
- lable in our papers on this subject
- Briefly, the function tries symmetry factorisation, fol-
- lowed by zero track elimination, followed by disconnect-
- ed subspace identifcation by path tracing.

## Implementation structure

- Symmetry and trajectory-level state space reduction. Tries all
- applicable reduction methods (unless disabled during the call
- to create.m) and returns a cell array of projectors into a set
- of independently evolving reduced subspaces. Syntax:
- projectors=reduce(spin_system,L,rho)
- L - Liouvillian matrix
- rho - initial state (source state screening) or
- destination state (destination state screening)
- projectors -a cell array of projectors into independently
- evolving reduced subspaces. The projectors are
- to be used as follows:
- L_reduced=P'*L*P; (for matrices)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `report()`, `iscell()`, `isfield()`, `true()`, `num2str()`, `irrep_keep_index()`, `zte()`, `path_trace()`.
