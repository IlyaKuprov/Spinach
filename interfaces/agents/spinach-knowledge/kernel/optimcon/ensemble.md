# kernel/optimcon/ensemble.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ensemble.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`
- Total lines: 295

## Purpose

A parallel wrapper around GRAPE that enables ensemble optimal control optimisations. This function handles systems with multiple control po- wer levels, multiple resonance offsets, multistate transfers, ensemb- les of drift Liouvillians, etc. Syntax: [traj_data,fidelity,... gradient,hessian]=ensemble(waveform,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `ens_block()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(spin_system,waveform)`.
- Lines 56-57: Worker-resident problem data handles; implemented by `invariants=spin_system.control.invariants`.
- Lines 60-62: Live problem data is the client-side control structure less what the workers already hold; implemented by `control=rmfield(spin_system.control,intersect({'invariants','drift_slices','worker_cases','basis'}, fieldnames(spin_system.control)))`.
- Lines 65-66: Count the outputs and the cases; implemented by `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 74-75: Run the ensemble loop, each worker over its own case block; implemented by `spmd (poolsize)`.
- Lines 77-79: Evaluate the block of cases assigned to this worker; implemented by `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value, control,spmdIndex,waveform,n_outputs)`.
- Lines 81-83: Reduce to the first worker and pack; implemented by `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1), 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 87-88: Collect from the first worker; implemented by `results=results{1}; traj_data=results.traj; fidelities=results.fid`.
- Lines 91-92: Average the block trajectory sums; implemented by `if ismember('average',control.traj_opts)`.
- Lines 100-101: Ensemble averages of fidelity, gradient, and Hessian; implemented by `fidelity=sum(fidelities)/n_cases`.
- Lines 109-110: Run diagnostic plotting (expensive!); implemented by `if ~isempty(spin_system.control.plotting)`.
- Lines 112-113: With or without instrumental distortions; implemented by `if ~isempty(spin_system.control.distplot)`.
- Lines 115-116: Apply the distortions; implemented by `dist_waveform=waveform`.
- Lines 119-120: Extract and apply distortion function; implemented by `dist_function=spin_system.control.distplot{k}`.
- Lines 125-126: Real-life trajectory and the distorted control sequence; implemented by `ctrl_trajan(spin_system,dist_waveform,traj_data,fidelities)`.
- Lines 130-131: Real-life trajectory but the ideal control sequence; implemented by `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`.

### Control flow inferred from the code

- Line 67: conditional branch on `(n_outputs>3)&&(~all(cellfun(@(f)isequal(f,@no_dist),control.distortion(:))))`.
- Line 70: conditional branch on `(n_outputs>3)&&(~strcmp(control.integrator,'rectangle'))`.
- Line 92: conditional branch on `ismember('average',control.traj_opts)`.
- Line 94: `for` loop over `n=2:numel(traj_data)`.
- Line 102: conditional branch on `n_outputs>2`.
- Line 105: conditional branch on `n_outputs>3`.
- Line 110: conditional branch on `~isempty(spin_system.control.plotting)`.
- Line 113: conditional branch on `~isempty(spin_system.control.distplot)`.
- Line 117: `for` loop over `k=1:numel(spin_system.control.distplot)`.

### Control flow inside `ens_block`

- Line 144: the live control structure is grafted over the frozen worker copy and every field the live copy lacks (frozen invariants, the case blocks, the waveform basis) is taken from the frozen copy; the drift generators come from the worker's own slice.
- Line 150: dispatch on `ss.bas.formalism` chooses the GRAPE function once, `grape_liouv` for `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, `grape_hilb` for `zeeman-hilb`.
- Line 169: `for` loop over the cases of the block, each case indexed through the six catalog columns (state pair, drift, power level, offset, phase-cycle step, distortion).
- Line 182: conditional branch on `~isempty(control.phase_cycle)`; the phase-cycle row multiplies the initial and target states by phase factors and rotates each control channel through the block-diagonal rotation matrix `R` built from the channel phases.
- Line 190: conditional branch on `~isempty(off_ens_sizes)`; `ind2sub` unpacks the offset combination index, first channel fastest, and `2*pi*offset` times each channel offset operator is added to the drift.
- Line 198: the waveform is scaled to physical units by the power level of the case.
- Line 202: `for` loop over the distortion functions, collecting their Jacobian only when derivatives are requested.
- Line 213: one GRAPE call with as many outputs as requested: trajectory, fidelity, gradient, Hessian.
- Line 218: the gradient passes through the distortion Jacobian and the transposed rotation, and is added to the block sum scaled by the power level.
- Line 224: the Hessian is rotated on both sides by the Kronecker product of the identity over time steps with the transposed rotation, and added to the block sum scaled by the squared power level.
- Line 231: conditional branch on `ismember('average',control.traj_opts)&&(n_mine>0)`; the forward trajectories of a non-empty block are added with the overloaded `plus`, which also covers the Hilbert-space cell trajectories, into one entry so that only block sums travel to the client, which divides by the case count; an empty block contributes nothing.

### Key state/data transformations

- Lines 57: computes `invariants` using `invariants=spin_system.control.invariants`.
- Lines 58: computes `drift_slices` using `drift_slices=spin_system.control.drift_slices`.
- Lines 61-62: computes `control` using `control=rmfield(spin_system.control,intersect({'invariants','drift_slices','worker_cases','basis'}, fieldnames(spin_system.control)))`.
- Lines 63: computes `control.return_traj` using `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 66: computes `n_outputs` using `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 78-79: computes `[traj_local,fid_local,grad_local,hess_local]` using `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value, control,spmdIndex,waveform,n_outputs)`.
- Lines 82-83: computes `results` using `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1), 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 89: computes `gradient` using `gradient=results.grad; hessian=results.hess`.
- Lines 93: computes `ave_traj` using `ave_traj=traj_data{1}.forward`.
- Lines 97: computes `traj_data` using `traj_data={struct('forward',{(1/n_cases)*ave_traj})}`.
- Lines 101: computes `fidelity` using `fidelity=sum(fidelities)/n_cases`.
- Lines 106: computes `hessian` using `hessian=reshape(hessian/n_cases,numel(waveform)*[1 1])`.
- Lines 116: computes `dist_waveform` using `dist_waveform=waveform`.
- Lines 120: computes `dist_function` using `dist_function=spin_system.control.distplot{k}`.

### Local helper functions

- Line 139: `ens_block()` — `function [traj,fid,grad,hess]=ens_block(ss,drifts,control,block,waveform,n_outputs)`. Graft live client data over the frozen worker copy, keep what only the worker holds
  - Representative operation: `frozen=ss.control; ss.control=control; ss.control.drifts=drifts`.
  - Representative operation: `missing=setdiff(fieldnames(frozen),fieldnames(control))`.
- Line 241: `grumble()` — `function grumble(spin_system,waveform)`.
  - Representative operation: `if ~isfield(spin_system,'control')`.
  - Representative operation: `error('control data missing from spin_system, run optimcon() first.')`.

## Parameters / inputs

- waveform -control coefficients for each control operator, rad/s

## Outputs

- traj_data -trajectory data for subsequent diagnostic plotting
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.
- gradient -gradient of the fidelity with respect to the control
- sequence. When penalty methods are specified, gradi-
- ent is returned as an array separating penalty gra-
- dients from the fidelity gradient.
- hessian -Hessian of the fidelity with respect to the control
- sequence. When penalty methods are specified, gradi-
- ent is returned as an array separating penalty Hes-
- sians from the fidelity Hessian.
- Note: the ensemble cases enumerated in spin_system.control.catalog
- are processed in the contiguous per-worker blocks assigned by
- optimcon.m in spin_system.control.worker_cases. Each worker
- holds the common frozen problem and the drift generators of its
- own block, published by optimcon.m as pool constants, and grafts
- the live client-side control structure on top of them, so only
- the waveform and the live control fields travel at each objec-
- tive evaluation; the gradient, the Hessian, and averaged trajec-
- tories are summed on the workers. This func-
- tion must be called from the client, on the pool that was
- open when optimcon.m ran: a worker holds only its own block.

## Implementation structure

- A parallel wrapper around GRAPE that enables ensemble optimal control
- optimisations. This function handles systems with multiple control po-
- wer levels, multiple resonance offsets, multistate transfers, ensemb-
- les of drift Liouvillians, etc. Syntax:
- [traj_data,fidelity,...
- gradient,hessian]=ensemble(waveform,spin_system)
- waveform -control coefficients for each control operator, rad/s
- traj_data -trajectory data for subsequent diagnostic plotting
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rmfield()`, `intersect()`, `fieldnames()`, `isfield()`, `cellfun()`, `isequal()`, `poolsize()`, `ens_block()`, `spmdCat()`, `spmdPlus()`, `ismember()`, `reshape()`, `strcmp()`, `ctrl_trajan()`, `setdiff()`, `kron()`, `ind2sub()`, `sparse()`, `speye()`, `grape_liouv()`, `grape_hilb()`, `getCurrentWorker()`, `gcp()`.
