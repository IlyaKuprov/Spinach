# kernel/reduce.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/reduce.m`
- Signature: `projectors=reduce(spin_system,L,rho)`
- Total lines: 315

## Purpose

Symmetry and trajectory-level state space reduction. Tries all applicable reduction methods (unless disabled during the call to create.m) and returns a cell array of projectors into a set of independently evolving reduced subspaces. Syntax: projectors=reduce(spin_system,L,rho)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

Horizontal stacks of wavefunctions or Liouville states retain their actual complex columns during symmetry screening. A sector is retained when any column exceeds the population threshold; Liouville stacks continue through zero-track elimination and path tracing without disabling reduction.

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
