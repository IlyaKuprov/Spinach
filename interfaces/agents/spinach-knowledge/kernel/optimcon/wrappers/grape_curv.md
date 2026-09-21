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
