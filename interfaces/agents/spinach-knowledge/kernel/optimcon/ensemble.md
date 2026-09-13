# kernel/optimcon/ensemble.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ensemble.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`
- Total lines: 504

## Purpose

A parallel wrapper around GRAPE that enables ensemble optimal control optimisations. This function handles systems with multiple control po- wer levels, multiple resonance offsets, multistate transfers, ensemb- les of drift Liouvillians, etc. Syntax: [traj_data,fidelity,... gradient,hessian]=ensemble(waveform,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `ens_block()`, `gcp()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 52-53: Check consistency; implemented by `grumble(spin_system,waveform)`.
- Lines 55-56: Pull the worker-resident problem data handle; implemented by `invariants=spin_system.control.invariants`.
- Lines 58-59: Live problem data is the client-side control structure; implemented by `control=rmfield(spin_system.control,'invariants')`.
- Lines 61-62: Default the trajectory return flag; implemented by `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 64-65: Count the outputs and the cases; implemented by `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 67-68: Run the ensemble loop, each worker over its own case block; implemented by `spmd (poolsize)`.
- Lines 70-73: Evaluate the block of cases assigned to this worker; implemented by `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,control, control.worker_cases{spmdIndex}, waveform,n_outputs)`.
- Lines 75-77: Reduce to the first worker and pack; implemented by `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',{spmdCat(fid_local,2,1)}, 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 81-82: Collect from the first worker; implemented by `results=results{1}; traj_data=results.traj; fidelities=results.fid`.
- Lines 85-86: Apply trajectory options; implemented by `if ismember('average',spin_system.control.traj_opts)`.
- Lines 88-89: Average the block trajectory sums; implemented by `ave_traj=(1/n_cases)*traj_data{1}.forward`.
- Lines 94-95: Overwrite traj_data; implemented by `traj_data=[]; traj_data{1}.forward=ave_traj`.
- Lines 99-100: Add up fidelities; implemented by `fidelities=cell2mat(fidelities)`.
- Lines 103-104: Normalise gradient; implemented by `if n_outputs>2`.
- Lines 108-109: Normalise Hessian; implemented by `if (n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')`.
- Lines 115-116: Run diagnostic plotting (expensive!); implemented by `if ~isempty(spin_system.control.plotting)`.
- Lines 118-119: With or without instrumental distortions; implemented by `if ~isempty(spin_system.control.distplot)`.
- Lines 121-122: Apply the distortions; implemented by `dist_waveform=waveform`.

### Control flow inferred from the code

- Line 86: conditional branch on `ismember('average',spin_system.control.traj_opts)`.
- Line 90: `for` loop over `n=2:numel(traj_data)`.
- Line 104: conditional branch on `n_outputs>2`.
- Line 109: conditional branch on `(n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')`.
- Line 116: conditional branch on `~isempty(spin_system.control.plotting)`.
- Line 119: conditional branch on `~isempty(spin_system.control.distplot)`.
- Line 123: `for` loop over `k=1:numel(spin_system.control.distplot)`.

### Control flow inside `ens_block`

- Line 170: `for` loop over the cases of the block, `m=1:n_mine`, each case indexed through the six catalog columns (state pair, drift, power level, offset, phase-cycle step, distortion).
- Line 185: conditional branch on `~isempty(control.phase_cycle)`; the initial state, the target state, and each complex control channel receive the phases of the current phase-cycle row.
- Line 220: conditional branch on `~isempty(off_ens_sizes)`; the linear offset index is unpacked into one index per offset channel and `2*pi*offset` times the channel offset operator is added to the drift.
- Line 245: the waveform is scaled to physical units by the power level of the case.
- Lines 248, 284, 335: dispatch on the number of outputs; two outputs apply the distortion functions and call GRAPE for fidelity and trajectory, three outputs also accumulate the distortion Jacobian and apply it to the gradient, four outputs call GRAPE for the Hessian as well.
- Line 262: dispatch on `ss.bas.formalism`; `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef` go to `grape_liouv`, `zeeman-hilb` goes to `grape_hilb`.
- Line 365: conditional branch on `(~isempty(control.phase_cycle))&&(n_outputs>2)`; the phase-cycle phases are removed from the gradient channel by channel.
- Line 389: conditional branch on `(~isempty(control.phase_cycle))&&(n_outputs>3)`; the Hessian is reshaped to `[ncont nsteps ncont nsteps]`, the phases are removed from both sides, and it is reshaped back.
- Line 429: the gradient and the Hessian of the case are scaled by the power level and its square and added to the running block sums.
- Line 440: conditional branch on `ismember('average',control.traj_opts)`; the forward trajectories of the block are summed into one entry so that only block sums travel to the client, which divides by the case count.

### Key state/data transformations

- Lines 56: computes `invariants` using `invariants=spin_system.control.invariants`.
- Lines 59: computes `control` using `control=rmfield(spin_system.control,'invariants')`.
- Lines 62: computes `control.return_traj` using `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 65: computes `n_outputs` using `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 71-73: computes `[traj_local,fid_local,grad_local,hess_local]` using `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,control, control.worker_cases{spmdIndex}, waveform,n_outputs)`.
- Lines 76-77: computes `results` using `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',{spmdCat(fid_local,2,1)}, 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 83: computes `gradient` using `gradient=results.grad; hessian=results.hess`.
- Lines 89: computes `ave_traj` using `ave_traj=(1/n_cases)*traj_data{1}.forward`.
- Lines 95: computes `traj_data` using `traj_data=[]; traj_data{1}.forward=ave_traj`.
- Lines 100: computes `fidelities` using `fidelities=cell2mat(fidelities)`.
- Lines 101: computes `fidelity` using `fidelity=sum(fidelities)/n_cases`.
- Lines 110: computes `hessian` using `hessian=reshape(hessian/n_cases,numel(waveform)*[1 1])`.
- Lines 122: computes `dist_waveform` using `dist_waveform=waveform`.
- Lines 126: computes `dist_function` using `dist_function=spin_system.control.distplot{k}`.

### Local helper functions

- Line 146: `ens_block()` — `function [traj,fid,grad,hess]=ens_block(ss,control,my_cases,waveform,n_outputs)`. Graft live client data over the frozen worker copy
  - Representative operation: `frozen=ss.control; ss.control=control`.
  - Representative operation: `for k=1:numel(control.frozen_fields)`.
- Line 451: `grumble()` — `function grumble(spin_system,waveform)`.
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
- holds the frozen problem data of its own block, published by
- optimcon.m as a pool constant, and grafts the live client-side
- control structure on top of it, so only the waveform and the
- live control fields travel at each objective evaluation; the
- gradient and the Hessian are summed on the workers. This func-
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

- Called routines detected from the main body: `grumble()`, `rmfield()`, `isfield()`, `poolsize()`, `ens_block()`, `spmdCat()`, `spmdPlus()`, `ismember()`, `cell2mat()`, `reshape()`, `strcmp()`, `dist_function()`, `ctrl_trajan()`, `cellfun()`, `fliplr()`, `cumprod()`, `sparse()`, `grape_liouv()`, `grape_hilb()`, `speye()`, `getCurrentWorker()`, `gcp()`.
