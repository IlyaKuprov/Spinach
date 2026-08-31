# kernel/optimcon/wrappers/grape_curv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/wrappers/grape_curv.m`
- Signature: `[traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,...`
- Total lines: 125

## Purpose

Cost function for optimal control using the GRAPE algorithm. Returns fidelity and gradient for a given waveform, specified in arbitrary curvilinear coordinates. Syntax: [traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,... dx_du,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Check consistency; implemented by `grumble(waveform_u,u2x,dx_du,spin_system)`.
- Lines 52-53: Count rectilinear and curvilinear coordinates; implemented by `size_x=spin_system.control.ncontrols`.
- Lines 57-58: Transform into rectilinear coordinates; implemented by `waveform_x=zeros(size_x,n_time_pts)`.
- Lines 63-64: Decide the output; implemented by `switch nargout`.
- Lines 66-67: Just fidelity; implemented by `case 2`.
- Lines 69-70: Call rectilinear coordinate GRAPE; implemented by `[traj_data,fidelity]=grape_xy(waveform_x,spin_system)`.
- Lines 72-73: Fidelity and gradient; implemented by `case 3`.
- Lines 75-76: Call rectilinear coordinate GRAPE; implemented by `[traj_data,fidelity,df_dx]=grape_xy(waveform_x,spin_system)`.
- Lines 78-79: Preallocate curvilinear gradients; implemented by `df_du=zeros(size_u,n_time_pts,size(df_dx,3))`.
- Lines 81-82: Fidelity and penalties; implemented by `for n=1:size(df_dx,3)`.
- Lines 84-85: Loop over time points; implemented by `for k=1:n_time_pts`.
- Lines 87-88: Translate gradients into curvilinear coordinates; implemented by `df_du(:,k,n)=dx_du(waveform_u(:,k))*df_dx(:,k,n)`.
- Lines 96-97: Complain and bomb out; implemented by `error('incorrect number of output arguments')`.

### Control flow inferred from the code

- Line 59: `for` loop over `n=1:n_time_pts`.
- Line 64: dispatches on `nargout`; cases `2`, `3`.
- Line 82: `for` loop over `n=1:size(df_dx,3)`.
- Line 85: `for` loop over `k=1:n_time_pts`.

### Key state/data transformations

- Lines 53: computes `size_x` using `size_x=spin_system.control.ncontrols`.
- Lines 54: computes `n_time_pts` using `n_time_pts=size(waveform_u,2)`.
- Lines 55: computes `size_u` using `size_u=size(waveform_u,1)`.
- Lines 58: computes `waveform_x` using `waveform_x=zeros(size_x,n_time_pts)`.
- Lines 60: computes `waveform_x(:,n)` using `waveform_x(:,n)=u2x(waveform_u(:,n))`.
- Lines 70: computes `[traj_data,fidelity]` using `[traj_data,fidelity]=grape_xy(waveform_x,spin_system)`.
- Lines 76: computes `[traj_data,fidelity,df_dx]` using `[traj_data,fidelity,df_dx]=grape_xy(waveform_x,spin_system)`.
- Lines 79: computes `df_du` using `df_du=zeros(size_u,n_time_pts,size(df_dx,3))`.
- Lines 88: computes `df_du(:,k,n)` using `df_du(:,k,n)=dx_du(waveform_u(:,k))*df_dx(:,k,n)`.

### Local helper functions

- Line 104: `grumble()` — `function grumble(waveform_u,u2x,dx_du,spin_system)`.
  - Representative operation: `if ~isfield(spin_system,'control')`.
  - Representative operation: `error('spin_system object lacks control information, run optimcon() first.')`.

## Parameters / inputs

- waveform_u -pulse waveform in curvilinear coordinates with indi-
- vidual coordinates in columns and time in rows
- u2x -a handle to a function that takes a column of curvi-
- linear coordinates and returns a column of coeffici-
- ents in front of the control operators
- dx_du -a handle to a function that takes a column of curvi-
- linear coordinates and returns the Jacobian matrix
- with the following structure:
- [dx(1)_du(1) dx(2)_du(1) dx(3)_du(1) ...
- dx(1)_du(2) dx(2)_du(2) dx(3)_du(2) ...
- ... ... ... ...]

## Outputs

- traj_data -system trajectory data structure used for visualisa-
- tion and progress reports
- fidelity -figure of merit for the overlap of the current state
- of the system and the desired state(s). When penalty
- methods are specified, fidelity is returned as an ar-
- ray separating the penalties from the simulation
- fidelity.
- df_du -gradient of the fidelity with respect to the control
- sequence. When penalty methods are specified, gradi-
- ent is returned as an array separating penalty gra-
- dients from the fidelity gradient.
- Note: penalities are computed using the rectilinear representation.

## Implementation structure

- Cost function for optimal control using the GRAPE algorithm. Returns
- fidelity and gradient for a given waveform, specified in arbitrary
- curvilinear coordinates. Syntax:
- [traj_data,fidelity,df_du]=grape_curv(waveform_u,u2x,...
- dx_du,spin_system)
- waveform_u - pulse waveform in curvilinear coordinates with indi-
- vidual coordinates in columns and time in rows
- u2x - a handle to a function that takes a column of curvi-
- linear coordinates and returns a column of coeffici-
- ents in front of the control operators
- dx_du - a handle to a function that takes a column of curvi-
- linear coordinates and returns the Jacobian matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `waveform_x()`, `u2x()`, `waveform_u()`, `grape_xy()`, `df_du()`, `dx_du()`, `df_dx()`, `isfield()`, `optimcon()`.
