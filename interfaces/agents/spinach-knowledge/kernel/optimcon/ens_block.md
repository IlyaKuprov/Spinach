# kernel/optimcon/ens_block.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_block.m`
- Signature: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`
- Total lines: 208

## Purpose

Fidelity, gradient, and Hessian contributions of one block of ensemble cases, evaluated on the parallel pool worker that holds the drift generators of that block. This function is called by ensemble.m inside its spmd block; the per-case physics (phase cycle, offsets, power level, waveform distortions, GRAPE) is applied here. Syntax: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

Freeze masks refer to the input waveform: the complete physical gradient is pulled back through all distortions, phase rotations, and power scaling before frozen input entries are zeroed. Supported Hessians likewise have zero frozen input rows and columns after the phase/power transformation; distortion Hessians remain unavailable. Direct GRAPE engine calls retain their separate existing behaviour.

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

## Control flow inside the case loop

- The live control structure is grafted over the frozen worker copy and every field the live copy lacks (frozen invariants, the case blocks, the waveform basis) is taken from the frozen copy; the drift generators come from the worker's own slice through the `drifts` argument.
- Dispatch on `spin_system.bas.formalism` chooses the GRAPE function once, `grape_liouv` for `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, `grape_hilb` for `zeeman-hilb`.
- `for` loop over the cases of the block, each case indexed through the six catalog columns (state pair, drift, power level, offset, phase-cycle step, distortion).
- Conditional branch on `~isempty(control.phase_cycle)`; the phase-cycle row multiplies the initial and target states by phase factors and rotates each control channel through the block-diagonal rotation matrix `R` built from the channel phases.
- Conditional branch on `~isempty(off_ens_sizes)`; `ind2sub` unpacks the offset combination index, first channel fastest, and `2*pi*offset` times each channel offset operator is added to the drift.
- The waveform is scaled to physical units by the power level of the case.
- `for` loop over the distortion functions, collecting their Jacobian only when derivatives are requested.
- One GRAPE call with as many outputs as requested: trajectory, fidelity, gradient, Hessian.
- The gradient passes through the distortion Jacobian and the transposed rotation, and is added to the block sum scaled by the power level.
- The Hessian is rotated on both sides by the Kronecker product of the identity over time steps with the transposed rotation, and added to the block sum scaled by the squared power level.
- Conditional branch on `ismember('average',control.traj_opts)&&(n_mine>0)`; the forward trajectories of a non-empty block are added with the overloaded `plus`, which also covers the Hilbert-space cell trajectories, into one entry so that only block sums travel to the client, which divides by the case count; an empty block contributes nothing.

## Header notes

The worker receives the frozen spin-system description and its local drift generators together with live control data, the block identifier, waveform, and requested output count.
