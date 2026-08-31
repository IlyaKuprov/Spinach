# kernel/optimcon/wrappers/grape_phase.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/wrappers/grape_phase.m`
- Signature: `[traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)`
- Total lines: 221

## Purpose

Cost function for optimal control using the GRAPE algorithm. Returns fidelity, gradient and Hessian for a given waveform, specified in po- lar coordinates. Only the phase channel gradient is returned, the am- plitude profile is taken as a given. Syntax: [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(spin_system,phi_profile)`.
- Lines 41-42: Extract the amplitude profile; implemented by `amplitudes=spin_system.control.amplitudes`.
- Lines 44-45: Count penalty terms; implemented by `npenterms=numel(spin_system.control.penalties)`.
- Lines 47-48: Transform into Cartesian coordinates; implemented by `waveform_xy=zeros(2*size(phi_profile,1),size(phi_profile,2))`.
- Lines 53-54: Translate the freeze mask to Cartesians; implemented by `spin_system.control.freeze=kron(spin_system.control.freeze,[1; 1])`.
- Lines 56-57: Just fidelity; implemented by `if nargout==2`.
- Lines 59-60: Call Cartesian GRAPE; implemented by `[traj_data,fidelity]=grape_xy(waveform_xy,spin_system)`.
- Lines 62-63: Fidelity and gradient; implemented by `elseif nargout==3`.
- Lines 65-66: Call Cartesian GRAPE; implemented by `[traj_data,fidelity,grad_xy]=grape_xy(waveform_xy,spin_system)`.
- Lines 68-69: Preallocate the answer; implemented by `gradient=zeros(size(phi_profile,1),size(phi_profile,2),npenterms+1)`.
- Lines 71-72: Loop over phase tracks; implemented by `for n=1:size(phi_profile,1)`.
- Lines 74-75: Loop over penalites; implemented by `for k=1:(npenterms+1)`.
- Lines 77-79: Translate derivatives; implemented by `[~,~,~,gradient(n,:,k)]=cartesian2polar(waveform_xy(2*n-1,:),waveform_xy(2*n,:), grad_xy(2*n-1,:,k), grad_xy(2*n,:,k))`.
- Lines 85-86: Fidelity, gradient and Hessian; implemented by `elseif nargout==4`.
- Lines 88-89: Call Cartesian ensemble GRAPE; implemented by `[traj_data,fidelity,grad_xy,hess_xy]=grape_xy(waveform_xy,spin_system)`.
- Lines 109-110: Transform Hessian from n^2 kxk block matrices to k^2 nxn matrices; implemented by `for n=1:size(hess_xy,3)`.
- Lines 119-120: Row-side waveform; implemented by `row_x=waveform_xy(2*n-1,:)`.
- Lines 123-124: Column-side waveform; implemented by `col_x=waveform_xy(2*k-1,:)`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:size(phi_profile,1)`.
- Line 57: conditional branch on `nargout==2`.
- Line 72: `for` loop over `n=1:size(phi_profile,1)`.
- Line 75: `for` loop over `k=1:(npenterms+1)`.
- Line 96: `for` loop over `n=1:size(phi_profile,1)`.
- Line 99: `for` loop over `k=1:(npenterms+1)`.
- Line 110: `for` loop over `n=1:size(hess_xy,3)`.
- Line 116: `for` loop over `n=1:size(phi_profile,1)`.
- Line 117: `for` loop over `k=1:size(phi_profile,1)`.
- Line 128: `for` loop over `m=1:(npenterms+1)`.
- Line 143: conditional branch on `n==k`.
- Line 158: `for` loop over `n=1:(npenterms+1)`.

### Key state/data transformations

- Lines 42: computes `amplitudes` using `amplitudes=spin_system.control.amplitudes`.
- Lines 45: computes `npenterms` using `npenterms=numel(spin_system.control.penalties)`.
- Lines 48: computes `waveform_xy` using `waveform_xy=zeros(2*size(phi_profile,1),size(phi_profile,2))`.
- Lines 50: computes `[waveform_xy(2*n-1,:),waveform_xy(2*n,:)]` using `[waveform_xy(2*n-1,:),waveform_xy(2*n,:)]=polar2cartesian(amplitudes(n,:),phi_profile(n,:))`.
- Lines 54: computes `spin_system.control.freeze` using `spin_system.control.freeze=kron(spin_system.control.freeze,[1; 1])`.
- Lines 60: computes `[traj_data,fidelity]` using `[traj_data,fidelity]=grape_xy(waveform_xy,spin_system)`.
- Lines 66: computes `[traj_data,fidelity,grad_xy]` using `[traj_data,fidelity,grad_xy]=grape_xy(waveform_xy,spin_system)`.
- Lines 69: computes `gradient` using `gradient=zeros(size(phi_profile,1),size(phi_profile,2),npenterms+1)`.
- Lines 78-79: computes `[~,~,~,gradient(n,:,k)]` using `[~,~,~,gradient(n,:,k)]=cartesian2polar(waveform_xy(2*n-1,:),waveform_xy(2*n,:), grad_xy(2*n-1,:,k), grad_xy(2*n,:,k))`.
- Lines 89: computes `[traj_data,fidelity,grad_xy,hess_xy]` using `[traj_data,fidelity,grad_xy,hess_xy]=grape_xy(waveform_xy,spin_system)`.
- Lines 93: computes `hessian` using `hessian=zeros(numel(phi_profile),numel(phi_profile),npenterms+1)`.
- Lines 111-112: computes `hess_xy(:,:,n)` using `hess_xy(:,:,n)=hess_reorder(hess_xy(:,:,n),spin_system.control.ncontrols, spin_system.control.pulse_nsteps)`.
- Lines 120: computes `row_x` using `row_x=waveform_xy(2*n-1,:)`.
- Lines 121: computes `row_y` using `row_y=waveform_xy(2*n,:)`.
- Lines 124: computes `col_x` using `col_x=waveform_xy(2*k-1,:)`.
- Lines 125: computes `col_y` using `col_y=waveform_xy(2*k,:)`.
- Lines 131: computes `Dx` using `Dx=grad_xy(2*n-1,:,m)`.
- Lines 132: computes `Dy` using `Dy=grad_xy(2*n,:,m)`.

### Local helper functions

- Line 168: `grumble()` — `function grumble(spin_system,phi_profile)`.
  - Representative operation: `if ~isfield(spin_system,'control')`.
  - Representative operation: `error('control data missing from spin_system, run optimcon() first.')`.

## Parameters / inputs

- phi_profile -set of control pulse phases from an amplitude-phase
- description.

## Outputs

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
- fidelity, gradient and Hessian for a given waveform, specified in po-
- lar coordinates. Only the phase channel gradient is returned, the am-
- plitude profile is taken as a given. Syntax:
- [traj_data,fidelity,gradient,hessian]=grape_phase(phi_profile,spin_system)
- phi_profile -set of control pulse phases from an amplitude-phase
- description.
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `waveform_xy()`, `polar2cartesian()`, `amplitudes()`, `phi_profile()`, `grape_xy()`, `gradient()`, `cartesian2polar()`, `grad_xy()`, `hess_xy()`, `hess_reorder()`, `hess_block()`, `hessian()`, `isfield()`, `optimcon()`, `any()`.
