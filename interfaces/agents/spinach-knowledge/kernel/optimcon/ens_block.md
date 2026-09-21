# kernel/optimcon/ens_block.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_block.m`
- Signature: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`
- Total lines: 198

## Purpose

Fidelity, gradient, and Hessian contributions of one block of ensemble cases, evaluated on the parallel pool worker that holds the drift generators of that block. This function is called by ensemble.m inside its spmd block; the per-case physics (phase cycle, offsets, power level, waveform distortions, GRAPE) is applied here. Syntax: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -frozen problem published by optimcon.m, with
- the drift generators removed
- drifts -cell array of drift generators, populated at
- the indices that the cases of this block use
- control -live client-side control structure
- block -index of the case block, into the cell array
- spin_system.control.worker_cases
- waveform -control coefficients for each control opera-
- tor, [ncontrols x nsteps], rad/s
- n_outputs -number of outputs requested from ensemble.m,
- 2 for the fidelity, 3 for the gradient, 4 for
- the Hessian

## Outputs

- traj -cell array of trajectory structures, one per
- case of the block in block order; when the
- control.traj_opts contains 'average', one
- structure holding the sum over the block, or
- an empty cell for an empty block
- fid -[1 x n_block] array of case fidelities
- grad -sum of the case gradients over the block, a
- [ncontrols*nsteps x 1] column, empty unless
- n_outputs>2
- hess -sum of the case Hessians over the block, a
- [(ncontrols*nsteps)^2 x 1] column, empty un-
- less n_outputs>3

## Implementation structure

- Fidelity, gradient, and Hessian contributions of one block of ensem-
- ble cases, evaluated on the parallel pool worker that holds the drift
- generators of that block. This function is called by ensemble.m in-
- side its spmd block; the per-case physics (phase cycle, offsets, po-
- wer level, waveform distortions, GRAPE) is applied here. Syntax:
- [traj,fid,grad,hess]=ens_block(spin_system,drifts,control,...
- block,waveform,n_outputs)
- spin_system -frozen problem published by optimcon.m, with
- the drift generators removed
- drifts -cell array of drift generators, populated at
- the indices that the cases of this block use
- control -live client-side control structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `setdiff()`, `fieldnames()`, `cellfun()`, `ind2sub()`, `sparse()`, `speye()`, `kron()`, `reshape()`, `grape_liouv()`, `grape_hilb()`, `ismember()`, `isfield()`.
