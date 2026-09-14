# kernel/optimcon/ens_block.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ens_block.m`
- Signature: `[traj,fid,grad,hess]=ens_block(spin_system,drifts,control,block,waveform,n_outputs)`
- Total lines: 186

## Purpose

Fidelity, gradient, and Hessian contributions of one block of ensem- ble cases, evaluated on the parallel pool worker that holds the drift generators of that block. This function is called by ensemble.m in- side its spmd block; the per-case physics (phase cycle, offsets, po- wer level, waveform distortions, GRAPE) is applied here. Syntax: [traj,fid,grad,hess]=ens_block(spin_system,drifts,control,... block,waveform,n_

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Check consistency; implemented by `grumble(spin_system,drifts,control,block,waveform,n_outputs)`.
- Lines 59-60: Graft live client data over the frozen worker copy, keep what only the worker holds; implemented by `frozen=spin_system.control; spin_system.control=control`.
- Lines 66-67: GRAPE function for the formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 76-77: Cases of this block, waveform dimensions, and offset ensemble size; implemented by `my_cases=frozen.worker_cases{block}; n_mine=numel(my_cases); catalog=control.catalog`.
- Lines 80-81: Preallocate block outputs, derivative buffers only when requested; implemented by `traj=cell(n_mine,1); fid=zeros(1,n_mine); grad=[]; hess=[]`.
- Lines 85-86: Loop over the cases of the block; implemented by `for m=1:n_mine`.
- Lines 88-89: Extract ensemble indices; implemented by `n=my_cases(m); n_rho=catalog(n,1); n_sys=catalog(n,2)`.
- Lines 93-94: Get initial and target states, drift, and waveform; implemented by `rho_init=control.rho_init{n_rho}; rho_targ=control.rho_targ{n_rho}`.
- Lines 97-98: Phase cycle: a rotation of each control channel and phases on the states; implemented by `R=eye(ncont)`.
- Lines 106-107: Add offset terms, first channel index fastest (user specifies offsets in Hz); implemented by `if ~isempty(off_ens_sizes)`.
- Lines 114-115: Move the waveform into physical units; implemented by `power_lvl=control.pwr_levels(n_pwr); local_waveform=power_lvl*local_waveform`.
- Lines 117-118: Apply waveform distortions, with their Jacobian when derivatives are needed; implemented by `if n_outputs>2, J=speye(numel(local_waveform)); end`.
- Lines 128-129: Fidelity, trajectory, and derivatives; implemented by `outputs=cell(1,n_outputs)`.
- Lines 133-134: Gradient through the Jacobian, the phase cycle, and the power level; implemented by `if n_outputs>2`.
- Lines 139-140: Hessian through the phase cycle and the power level; implemented by `if n_outputs>3`.
- Lines 147-148: Collapse a non-empty block into one trajectory sum when only the average is needed; implemented by `if ismember('average',control.traj_opts)&&(n_mine>0)`.

### Control flow inferred from the code

- Line 62: `for` loop over `k=1:numel(missing)`.
- Line 67: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv','zeeman-wavef'}`, `'zeeman-hilb'`.
- Line 82: conditional branch on `n_outputs>2, grad=zeros(ncont*nsteps,1); end`.
- Line 83: conditional branch on `n_outputs>3, hess=zeros((ncont*nsteps)^2,1); end`.
- Line 86: `for` loop over `m=1:n_mine`.
- Line 99: conditional branch on `~isempty(control.phase_cycle)`.
- Line 107: conditional branch on `~isempty(off_ens_sizes)`.
- Line 109: `for` loop over `k=1:numel(off_ens_sizes)`.
- Line 118: conditional branch on `n_outputs>2, J=speye(numel(local_waveform)); end`.
- Line 119: `for` loop over `k=1:size(control.distortion,2)`.
- Line 120: conditional branch on `n_outputs>2`.
- Line 134: conditional branch on `n_outputs>2`.
- Line 140: conditional branch on `n_outputs>3`.
- Line 148: conditional branch on `ismember('average',control.traj_opts)&&(n_mine>0)`.

### Control flow inside the case loop

- Line 61: the live control structure is grafted over the frozen worker copy and every field the live copy lacks (frozen invariants, the case blocks, the waveform basis) is taken from the frozen copy; the drift generators come from the worker's own slice through the `drifts` argument.
- Line 67: dispatch on `spin_system.bas.formalism` chooses the GRAPE function once, `grape_liouv` for `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, `grape_hilb` for `zeeman-hilb`.
- Line 86: `for` loop over the cases of the block, each case indexed through the six catalog columns (state pair, drift, power level, offset, phase-cycle step, distortion).
- Line 99: conditional branch on `~isempty(control.phase_cycle)`; the phase-cycle row multiplies the initial and target states by phase factors and rotates each control channel through the block-diagonal rotation matrix `R` built from the channel phases.
- Line 107: conditional branch on `~isempty(off_ens_sizes)`; `ind2sub` unpacks the offset combination index, first channel fastest, and `2*pi*offset` times each channel offset operator is added to the drift.
- Line 115: the waveform is scaled to physical units by the power level of the case.
- Line 119: `for` loop over the distortion functions, collecting their Jacobian only when derivatives are requested.
- Line 130: one GRAPE call with as many outputs as requested: trajectory, fidelity, gradient, Hessian.
- Line 135: the gradient passes through the distortion Jacobian and the transposed rotation, and is added to the block sum scaled by the power level.
- Line 141: the Hessian is rotated on both sides by the Kronecker product of the identity over time steps with the transposed rotation, and added to the block sum scaled by the squared power level.
- Line 148: conditional branch on `ismember('average',control.traj_opts)&&(n_mine>0)`; the forward trajectories of a non-empty block are added with the overloaded `plus`, which also covers the Hilbert-space cell trajectories, into one entry so that only block sums travel to the client, which divides by the case count; an empty block contributes nothing.

### Key state/data transformations

- Lines 60: computes `frozen` using `frozen=spin_system.control; spin_system.control=control`.
- Lines 61: computes `missing` using `missing=setdiff(fieldnames(frozen),fieldnames(control))`.
- Lines 63: computes `spin_system.control.(missing{k})` using `spin_system.control.(missing{k})=frozen.(missing{k})`.
- Lines 69: computes `grape` using `grape=@grape_liouv`.
- Lines 77: computes `my_cases` using `my_cases=frozen.worker_cases{block}; n_mine=numel(my_cases); catalog=control.catalog`.
- Lines 78: computes `ncont` using `ncont=size(waveform,1); nsteps=size(waveform,2); off_ens_sizes=cellfun(@numel,control.offsets)`.
- Lines 81: computes `traj` using `traj=cell(n_mine,1); fid=zeros(1,n_mine); grad=[]; hess=[]`.
- Lines 89: computes `n` using `n=my_cases(m); n_rho=catalog(n,1); n_sys=catalog(n,2)`.
- Lines 90: computes `n_pwr` using `n_pwr=catalog(n,3); n_off=catalog(n,4)`.
- Lines 91: computes `n_phi` using `n_phi=catalog(n,5); n_dis=catalog(n,6)`.
- Lines 94: computes `rho_init` using `rho_init=control.rho_init{n_rho}; rho_targ=control.rho_targ{n_rho}`.
- Lines 95: computes `L` using `L=drifts{n_sys}; local_waveform=waveform`.
- Lines 98: computes `R` using `R=eye(ncont)`.
- Lines 100: computes `phi` using `phi=control.phase_cycle(n_phi,:)`.
- Lines 103: computes `local_waveform` using `local_waveform=R*local_waveform`.
- Lines 108: computes `off_idx` using `off_idx=cell(1,numel(off_ens_sizes)); [off_idx{:}]=ind2sub([off_ens_sizes 1],n_off)`.
- Lines 115: computes `power_lvl` using `power_lvl=control.pwr_levels(n_pwr); local_waveform=power_lvl*local_waveform`.
- Lines 121: computes `[local_waveform,stage_jacobian]` using `[local_waveform,stage_jacobian]=control.distortion{n_dis,k}(local_waveform)`.

### Local helper functions

- Line 159: `grumble()` — `function grumble(spin_system,drifts,control,block,waveform,n_outputs)`.
  - Representative operation: `if (~isfield(spin_system,'control'))||(~isfield(spin_system.control,'worker_cases'))`.
  - Representative operation: `error('spin_system must be the frozen problem published by optimcon().')`.

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
- 1 for the trajectory, 2 for the fidelity, 3
- for the gradient, 4 for the Hessian

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
