# kernel/optimcon/ensemble.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ensemble.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`
- Total lines: 488

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

- Lines 50-51: Check consistency; implemented by `grumble(spin_system,waveform)`.
- Lines 53-54: Pull the worker-resident problem data handle; implemented by `invariants=spin_system.control.invariants`.
- Lines 56-57: Live problem data is the client-side control structure; implemented by `control=rmfield(spin_system.control,'invariants')`.
- Lines 59-60: Default the trajectory return flag; implemented by `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 62-63: Count the outputs and the cases; implemented by `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 65-66: Run the ensemble loop, each worker over its own case block; implemented by `spmd (poolsize)`.
- Lines 68-71: Evaluate the block of cases assigned to this worker; implemented by `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,control, control.worker_cases{spmdIndex}, waveform,n_outputs)`.
- Lines 73-75: Reduce to the first worker and pack; implemented by `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',{spmdCat(fid_local,2,1)}, 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 79-80: Collect from the first worker; implemented by `results=results{1}; traj_data=results.traj; fidelities=results.fid`.
- Lines 83-84: Apply trajectory options; implemented by `if ismember('average',spin_system.control.traj_opts)`.
- Lines 86-87: Return average trajectory; implemented by `ave_traj=(1/n_cases)*traj_data{1}.forward`.
- Lines 92-93: Overwrite traj_data; implemented by `traj_data=[]; traj_data{1}.forward=ave_traj`.
- Lines 97-98: Add up fidelities; implemented by `fidelities=cell2mat(fidelities)`.
- Lines 101-102: Normalise gradient; implemented by `if n_outputs>2`.
- Lines 106-107: Normalise Hessian; implemented by `if (n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')`.
- Lines 113-114: Run diagnostic plotting (expensive!); implemented by `if ~isempty(spin_system.control.plotting)`.
- Lines 116-117: With or without instrumental distortions; implemented by `if ~isempty(spin_system.control.distplot)`.
- Lines 119-120: Apply the distortions; implemented by `dist_waveform=waveform`.

### Control flow inferred from the code

- Line 84: conditional branch on `ismember('average',spin_system.control.traj_opts)`.
- Line 88: `for` loop over `n=2:numel(traj_data)`.
- Line 102: conditional branch on `n_outputs>2`.
- Line 107: conditional branch on `(n_outputs>3)&&strcmp(spin_system.control.integrator,'rectangle')`.
- Line 114: conditional branch on `~isempty(spin_system.control.plotting)`.
- Line 117: conditional branch on `~isempty(spin_system.control.distplot)`.
- Line 121: `for` loop over `k=1:numel(spin_system.control.distplot)`.

### Key state/data transformations

- Lines 54: computes `invariants` using `invariants=spin_system.control.invariants`.
- Lines 57: computes `control` using `control=rmfield(spin_system.control,'invariants')`.
- Lines 60: computes `control.return_traj` using `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 63: computes `n_outputs` using `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 69-71: computes `[traj_local,fid_local,grad_local,hess_local]` using `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,control, control.worker_cases{spmdIndex}, waveform,n_outputs)`.
- Lines 74-75: computes `results` using `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',{spmdCat(fid_local,2,1)}, 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 81: computes `gradient` using `gradient=results.grad; hessian=results.hess`.
- Lines 87: computes `ave_traj` using `ave_traj=(1/n_cases)*traj_data{1}.forward`.
- Lines 93: computes `traj_data` using `traj_data=[]; traj_data{1}.forward=ave_traj`.
- Lines 98: computes `fidelities` using `fidelities=cell2mat(fidelities)`.
- Lines 99: computes `fidelity` using `fidelity=sum(fidelities)/n_cases`.
- Lines 108: computes `hessian` using `hessian=reshape(hessian/n_cases,numel(waveform)*[1 1])`.
- Lines 120: computes `dist_waveform` using `dist_waveform=waveform`.
- Lines 124: computes `dist_function` using `dist_function=spin_system.control.distplot{k}`.

### Local helper functions

- Line 144: `ens_block()` — `function [traj,fid,grad,hess]=ens_block(ss,control,my_cases,waveform,n_outputs)`. Graft live client data over the frozen worker copy
  - Representative operation: `frozen=ss.control; ss.control=control`.
  - Representative operation: `for k=1:numel(control.frozen_fields)`.
- Line 440: `grumble()` — `function grumble(spin_system,waveform)`.
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
- gradient and the Hessian are summed on the workers.

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

- Called routines detected from the main body: `grumble()`, `rmfield()`, `isfield()`, `poolsize()`, `ens_block()`, `spmdCat()`, `spmdPlus()`, `ismember()`, `cell2mat()`, `reshape()`, `strcmp()`, `dist_function()`, `ctrl_trajan()`, `cellfun()`, `fliplr()`, `cumprod()`, `sparse()`, `grape_liouv()`, `grape_hilb()`, `speye()`.
