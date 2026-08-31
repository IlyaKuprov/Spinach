# kernel/optimcon/wrappers/grape_xy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/wrappers/grape_xy.m`
- Signature: `[traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)`
- Total lines: 213

## Purpose

Cost function for optimal control using the GRAPE algorithm. Returns fidelity, gradient and hessian for a given waveform, specified in Car- tesian coordinates (x and y channels). Syntax: [traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,waveform)`.
- Lines 41-42: Count penalty terms; implemented by `npenterms=numel(spin_system.control.penalties)`.
- Lines 44-46: Freeze masks act on time points, waveform bases mix them; implemented by `if (~isempty(spin_system.control.basis))&& (~isempty(spin_system.control.freeze))`.
- Lines 50-51: Translate the basis if necessary; implemented by `if ~isempty(spin_system.control.basis)`.
- Lines 56-57: Decide how much needs computing; implemented by `if nargout==2`.
- Lines 59-60: Preallocate output; implemented by `fidelity=zeros(1,npenterms+1)`.
- Lines 62-63: Calculate the objective and its gradient; implemented by `[traj_data,fidelity(1)]=ensemble(waveform,spin_system)`.
- Lines 65-66: Apply all penalties; implemented by `for n=1:npenterms`.
- Lines 68-71: Calculate the penalty (this should be moved to ensemble); implemented by `pen=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 73-74: Add to fidelity array; implemented by `fidelity(n+1)=spin_system.control.p_weights(n)*pen`.
- Lines 84-85: Calculate the objective and its gradient; implemented by `[traj_data,fidelity(1),grad(:,:,1)]=ensemble(waveform,spin_system)`.
- Lines 90-93: Calculate the penalty; implemented by `[pen,pen_grad]=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 95-96: Add to relevant arrays; implemented by `fidelity(n+1)=spin_system.control.p_weights(n)*pen`.
- Lines 114-115: Calculate objective, gradient and Hessian; implemented by `[traj_data,fidelity(1),grad(:,:,1),hess(:,:,1)]=ensemble(waveform,spin_system)`.
- Lines 120-123: Calculate the penalty; implemented by `[pen,pen_grad,pen_hess]=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 125-126: Add to relevant arrays; implemented by `fidelity(:,n+1)=spin_system.control.p_weights(n)*pen`.
- Lines 135-136: Transform gradients; implemented by `grad=tensorprod(spin_system.control.basis,grad,2,2)`.
- Lines 139-142: Preallocate transformed Hessians; implemented by `hess_in_basis=zeros([nwaves*spin_system.control.ncontrols nwaves*spin_system.control.ncontrols npenterms+1])`.

### Control flow inferred from the code

- Line 45: conditional branch on `(~isempty(spin_system.control.basis))&&`.
- Line 51: conditional branch on `~isempty(spin_system.control.basis)`.
- Line 57: conditional branch on `nargout==2`.
- Line 66: `for` loop over `n=1:npenterms`.
- Line 88: `for` loop over `n=1:npenterms`.
- Line 102: conditional branch on `~isempty(spin_system.control.basis)`.
- Line 118: `for` loop over `n=1:npenterms`.
- Line 133: conditional branch on `~isempty(spin_system.control.basis)`.
- Line 145: `for` loop over `n=1:size(hess,3)`.

### Key state/data transformations

- Lines 42: computes `npenterms` using `npenterms=numel(spin_system.control.penalties)`.
- Lines 52: computes `nwaves` using `nwaves=size(spin_system.control.basis,1)`.
- Lines 53: computes `waveform` using `waveform=waveform*spin_system.control.basis`.
- Lines 60: computes `fidelity` using `fidelity=zeros(1,npenterms+1)`.
- Lines 63: computes `[traj_data,fidelity(1)]` using `[traj_data,fidelity(1)]=ensemble(waveform,spin_system)`.
- Lines 69-71: computes `pen` using `pen=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 74: computes `fidelity(n+1)` using `fidelity(n+1)=spin_system.control.p_weights(n)*pen`.
- Lines 82: computes `grad` using `grad=zeros(size(waveform,1),size(waveform,2),npenterms+1)`.
- Lines 85: computes `[traj_data,fidelity(1),grad(:,:,1)]` using `[traj_data,fidelity(1),grad(:,:,1)]=ensemble(waveform,spin_system)`.
- Lines 91-93: computes `[pen,pen_grad]` using `[pen,pen_grad]=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 97: computes `grad(:,:,n+1)` using `grad(:,:,n+1)=spin_system.control.p_weights(n)*pen_grad`.
- Lines 112: computes `hess` using `hess=zeros(numel(waveform),numel(waveform),npenterms+1)`.
- Lines 115: computes `[traj_data,fidelity(1),grad(:,:,1),hess(:,:,1)]` using `[traj_data,fidelity(1),grad(:,:,1),hess(:,:,1)]=ensemble(waveform,spin_system)`.
- Lines 121-123: computes `[pen,pen_grad,pen_hess]` using `[pen,pen_grad,pen_hess]=penalty(waveform,spin_system.control.penalties{n}, spin_system.control.l_bound, spin_system.control.u_bound)`.
- Lines 126: computes `fidelity(:,n+1)` using `fidelity(:,n+1)=spin_system.control.p_weights(n)*pen`.
- Lines 128: computes `hess(:,:,n+1)` using `hess(:,:,n+1)=spin_system.control.p_weights(n)*pen_hess`.
- Lines 140-142: computes `hess_in_basis` using `hess_in_basis=zeros([nwaves*spin_system.control.ncontrols nwaves*spin_system.control.ncontrols npenterms+1])`.
- Lines 148-149: computes `hess_re` using `hess_re=hess_reorder(hess(:,:,n),spin_system.control.ncontrols, spin_system.control.pulse_nsteps)`.

### Local helper functions

- Line 182: `grumble()` — `function grumble(spin_system,waveform)`.
  - Representative operation: `if ~isfield(spin_system,'control')`.
  - Representative operation: `error('control data missing from spin_system, run optimcon() first.')`.

## Parameters / inputs

- waveform -normalised set of control amplitudes.

## Outputs

- traj_data -trajectory data
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

## Implementation structure

- Cost function for optimal control using the GRAPE algorithm. Returns
- fidelity, gradient and hessian for a given waveform, specified in Car-
- tesian coordinates (x and y channels). Syntax:
- [traj_data,fidelity,grad,hess]=grape_xy(waveform,spin_system)
- waveform -normalised set of control amplitudes.
- traj_data -trajectory data
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.
- gradient -gradient of the fidelity with respect to the control

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fidelity()`, `ensemble()`, `penalty()`, `grad()`, `tensorprod()`, `hess()`, `hess_reorder()`, `hess_in_basis()`, `isfield()`, `optimcon()`, `int2str()`.
