# kernel/optimcon/ensemble.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ensemble.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=ensemble(waveform,spin_system)`
- Total lines: 468

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

- Lines 48-49: Check consistency; implemented by `grumble(spin_system,waveform)`.
- Lines 51-52: Pull the ensemble case catalog; implemented by `catalog=spin_system.control.catalog`.
- Lines 55-56: Pull the worker-resident problem data handle; implemented by `invariants=spin_system.control.invariants`.
- Lines 58-59: Live problem data is the client-side control structure; implemented by `control=rmfield(spin_system.control,'invariants')`.
- Lines 61-62: Default the trajectory return flag; implemented by `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 64-65: Get offset ensemble size; implemented by `off_ens_sizes=cellfun(@numel,control.offsets)`.
- Lines 67-68: Count the outputs; implemented by `n_outputs=nargout`.
- Lines 70-71: Preallocate outputs; implemented by `traj_data=cell(n_cases,1); fidelities=cell(1,n_cases)`.
- Lines 74-75: Waveform dimension statistics; implemented by `ncont=size(waveform,1); nsteps=size(waveform,2)`.
- Lines 77-78: Parallelise over the ensemble; implemented by `nworkers=poolsize`.
- Lines 80-81: Run the ensemble loop; implemented by `parfor (n=1:n_cases,nworkers)`.
- Lines 83-84: Fetch worker-resident problem data; implemented by `ss=invariants.Value`.
- Lines 86-87: Graft live client data over the frozen worker copy; implemented by `frozen=ss.control; ss.control=control`.
- Lines 93-94: Extract ensemble indices; implemented by `n_rho=catalog(n,1); n_sys=catalog(n,2)`.
- Lines 98-99: Get initial and target state; implemented by `rho_init=control.rho_init{n_rho}`.
- Lines 102-103: Localise the waveform; implemented by `local_waveform=waveform`.
- Lines 105-106: Apply the phase cycle; implemented by `if ~isempty(control.phase_cycle)`.
- Lines 108-109: Apply phase to the initial state; implemented by `phi=control.phase_cycle(n_phi,1)`.

### Control flow inferred from the code

- Line 81: `parfor` loop over `(n=1:n_cases,nworkers)`.
- Line 88: `for` loop over `k=1:numel(control.frozen_fields)`.
- Line 106: conditional branch on `~isempty(control.phase_cycle)`.
- Line 117: `for` loop over `k=1:(size(local_waveform,1)/2)`.
- Line 141: conditional branch on `~isempty(off_ens_sizes)`.
- Line 146: `for` loop over `k=1:numel(off_ens_sizes)`.
- Line 169: conditional branch on `n_outputs==2`.
- Line 172: `for` loop over `k=1:size(control.distortion,2)`.
- Line 183: dispatches on `ss.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv','zeeman-wavef'}`, `'zeeman-hilb'`.
- Line 211: `for` loop over `k=1:size(control.distortion,2)`.
- Line 225: dispatches on `ss.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv','zeeman-wavef'}`, `'zeeman-hilb'`.
- Line 259: dispatches on `ss.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv','zeeman-wavef'}`, `'zeeman-hilb'`.
- Line 286: conditional branch on `(~isempty(control.phase_cycle))&&(n_outputs>2)`.
- Line 289: `for` loop over `k=1:(size(gradients{n},1)/2)`.

### Key state/data transformations

- Lines 52: computes `catalog` using `catalog=spin_system.control.catalog`.
- Lines 53: computes `n_cases` using `n_cases=size(catalog,1)`.
- Lines 56: computes `invariants` using `invariants=spin_system.control.invariants`.
- Lines 59: computes `control` using `control=rmfield(spin_system.control,'invariants')`.
- Lines 62: computes `control.return_traj` using `control.return_traj=isfield(control,'return_traj')&&control.return_traj`.
- Lines 65: computes `off_ens_sizes` using `off_ens_sizes=cellfun(@numel,control.offsets)`.
- Lines 68: computes `n_outputs` using `n_outputs=nargout`.
- Lines 71: computes `traj_data` using `traj_data=cell(n_cases,1); fidelities=cell(1,n_cases)`.
- Lines 72: computes `gradients` using `gradients=cell(1,n_cases); hessians=cell(1,n_cases)`.
- Lines 75: computes `ncont` using `ncont=size(waveform,1); nsteps=size(waveform,2)`.
- Lines 78: computes `nworkers` using `nworkers=poolsize`.
- Lines 84: computes `ss` using `ss=invariants.Value`.
- Lines 87: computes `frozen` using `frozen=ss.control; ss.control=control`.
- Lines 89: computes `fname` using `fname=control.frozen_fields{k}`.
- Lines 90: computes `ss.control.(fname)` using `ss.control.(fname)=frozen.(fname)`.
- Lines 94: computes `n_rho` using `n_rho=catalog(n,1); n_sys=catalog(n,2)`.
- Lines 95: computes `n_pwr` using `n_pwr=catalog(n,3); n_off=catalog(n,4)`.
- Lines 96: computes `n_phi` using `n_phi=catalog(n,5); n_dis=catalog(n,6)`.

### Local helper functions

- Line 423: `grumble()` — `function grumble(spin_system,waveform)`.
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
- are distributed over the parallel pool workers. Each case
- fetches the frozen problem data from the pool constant pub-
- lished by optimcon.m and grafts the live client-side control
- structure on top of it, so only the waveform and the live
- control fields travel at each objective evaluation.

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

- Called routines detected from the main body: `grumble()`, `rmfield()`, `isfield()`, `cellfun()`, `parfor()`, `catalog()`, `local_waveform()`, `fliplr()`, `cumprod()`, `cum_sizes()`, `dist_function()`, `grape_liouv()`, `grape_hilb()`, `speye()`, `ismember()`, `cell2mat()`.
