# kernel/optimcon/ensemble.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ensemble.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`
- Total lines: 192

## Purpose

A parallel wrapper around GRAPE that enables ensemble optimal control optimisations. This function handles systems with multiple control po- wer levels, multiple resonance offsets, multistate transfers, ensemb- les of drift Liouvillians, etc. Syntax: [traj_data,fidelity,... gradient,hessian]=ensemble(waveform,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Check consistency; implemented by `grumble(spin_system,waveform,nargout)`.
- Lines 56-57: Worker-resident problem data handles; implemented by `invariants=spin_system.control.invariants`.
- Lines 60-61: Live problem data is the client-side control structure less what the workers already hold; implemented by `control=rmfield(spin_system.control,{'invariants','drift_slices','worker_cases','basis'})`.
- Lines 64-65: Count the outputs and the cases; implemented by `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 67-68: Run the ensemble loop, each worker over its own case block; implemented by `spmd (poolsize)`.
- Lines 70-72: Evaluate the block of cases assigned to this worker; implemented by `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value, control,spmdIndex,waveform,n_outputs)`.
- Lines 74-76: Reduce to the first worker and pack; implemented by `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1), 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 80-81: Collect from the first worker, fidelities back into catalog order; implemented by `results=results{1}; gradient=results.grad; hessian=results.hess`.
- Lines 85-86: Average the block trajectory sums, or put the trajectories back into catalog order; implemented by `if ismember('average',control.traj_opts)`.
- Lines 96-97: Ensemble averages of fidelity, gradient, and Hessian; implemented by `fidelity=sum(fidelities)/n_cases`.
- Lines 105-106: Run diagnostic plotting (expensive!); implemented by `if ~isempty(spin_system.control.plotting)`.
- Lines 108-109: With or without instrumental distortions; implemented by `if ~isempty(spin_system.control.distplot)`.
- Lines 111-112: Apply the distortions; implemented by `dist_waveform=waveform`.
- Lines 115-116: Extract and apply distortion function; implemented by `dist_function=spin_system.control.distplot{k}`.
- Lines 121-122: Real-life trajectory and the distorted control sequence; implemented by `ctrl_trajan(spin_system,dist_waveform,traj_data,fidelities)`.
- Lines 126-127: Real-life trajectory but the ideal control sequence; implemented by `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`.

### Control flow inferred from the code

- Line 86: conditional branch on `ismember('average',control.traj_opts)`.
- Line 88: `for` loop over `n=2:numel(results.traj)`.
- Line 98: conditional branch on `n_outputs>2`.
- Line 101: conditional branch on `n_outputs>3`.
- Line 106: conditional branch on `~isempty(spin_system.control.plotting)`.
- Line 109: conditional branch on `~isempty(spin_system.control.distplot)`.
- Line 113: `for` loop over `k=1:numel(spin_system.control.distplot)`.

### Key state/data transformations

- Lines 57: computes `invariants` using `invariants=spin_system.control.invariants`.
- Lines 58: computes `drift_slices` using `drift_slices=spin_system.control.drift_slices`.
- Lines 61: computes `control` using `control=rmfield(spin_system.control,{'invariants','drift_slices','worker_cases','basis'})`.
- Lines 62: computes `control.return_traj` using `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 65: computes `n_outputs` using `n_outputs=nargout; n_cases=size(control.catalog,1)`.
- Lines 71-72: computes `[traj_local,fid_local,grad_local,hess_local]` using `[traj_local,fid_local,grad_local,hess_local]=ens_block(invariants.Value,drift_slices.Value, control,spmdIndex,waveform,n_outputs)`.
- Lines 75-76: computes `results` using `results=struct('traj',{spmdCat(traj_local,1,1)},'fid',spmdCat(fid_local,2,1), 'grad',spmdPlus(grad_local,1),'hess',spmdPlus(hess_local,1))`.
- Lines 82: computes `order` using `order=[spin_system.control.worker_cases{:}]`.
- Lines 83: computes `fidelities` using `fidelities=zeros(1,n_cases); fidelities(order)=results.fid`.
- Lines 87: computes `ave_traj` using `ave_traj=results.traj{1}.forward`.
- Lines 91: computes `traj_data` using `traj_data={struct('forward',{(1/n_cases)*ave_traj})}`.
- Lines 97: computes `fidelity` using `fidelity=sum(fidelities)/n_cases`.
- Lines 99: computes `gradient` using `gradient=reshape(gradient/n_cases,size(waveform))`.
- Lines 102: computes `hessian` using `hessian=reshape(hessian/n_cases,numel(waveform)*[1 1])`.
- Lines 112: computes `dist_waveform` using `dist_waveform=waveform`.
- Lines 116: computes `dist_function` using `dist_function=spin_system.control.distplot{k}`.

### Local helper functions

- Line 136: `grumble()` — `function grumble(spin_system,waveform,n_outputs)`.
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
- own block, published by optimcon.m as pool constants; the per-
- case physics runs in ens_block.m on each worker, so only the
- waveform and the live control fields travel at each objective
- evaluation, and the gradient, the Hessian, and averaged trajec-
- tories are summed on the workers. This function must be called
- from the client, on the pool that was open when optimcon.m ran:
- a worker holds only its own block.

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

- Called routines detected from the main body: `grumble()`, `rmfield()`, `isfield()`, `poolsize()`, `ens_block()`, `spmdCat()`, `spmdPlus()`, `ismember()`, `reshape()`, `ctrl_trajan()`, `gcp()`, `getCurrentWorker()`, `cellfun()`, `isequal()`, `strcmp()`, `ens_catalog()`.
