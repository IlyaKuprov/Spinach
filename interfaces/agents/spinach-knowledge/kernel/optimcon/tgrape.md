# kernel/optimcon/tgrape.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/tgrape.m`
- Signature: `[fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...`
- Total lines: 168

## Purpose

A special case of Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient with respect to the vector of waveform slice du- rations. Syntax: [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,... dt_grid,time_unit,rho_init,rho_targ)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- drift -the drift Liouvillian, a matrix
- controls -control operators, a cell array
- of matrices
- waveform -control coefficients for each control
- operator (columns) at each time slice
- (rows), rad/s
- dt_grid -time slice durations, a col vector
- in the units of time chosen so that
- the elements are of the order of 1
- time_unit -unit of time, seconds; this is needed
- because optimisers get stuck when the
- variables are badly scaled
- rho_init -initial state of the system, a column
- vector
- rho_targ -target state of the system, a column
- vector

## Outputs

- fidelity -fidelity of the control sequence
- grad -gradient of the fidelity with respect to
- the durations of the slices

## Implementation structure

- A special case of Gradient Ascent Pulse Engineering (GRAPE) objective
- function and gradient with respect to the vector of waveform slice du-
- rations. Syntax:
- [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...
- dt_grid,time_unit,rho_init,rho_targ)
- drift -the drift Liouvillian, a matrix
- controls -control operators, a cell array
- of matrices
- waveform -control coefficients for each control
- operator (columns) at each time slice
- (rows), rad/s
- dt_grid -time slice durations, a col vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fwd_traj()`, `bwd_traj()`, `waveform()`, `step()`, `dt_grid()`, `fliplr()`, `grad()`, `ismember()`, `iscolumn()`, `iscell()`, `any()`.
