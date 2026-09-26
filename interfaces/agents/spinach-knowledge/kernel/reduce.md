# kernel/reduce.m

- Signature: `projectors=reduce(spin_system,L,rho)`

## Purpose

Symmetry and trajectory-level state space reduction. Tries all applicable reduction methods (unless disabled during the call to create.m) and returns a cell array of projectors into a set of independently evolving reduced subspaces. Syntax: projectors=reduce(spin_system,L,rho)

## Physical / mathematical content

## Numerical / algorithmic content

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
