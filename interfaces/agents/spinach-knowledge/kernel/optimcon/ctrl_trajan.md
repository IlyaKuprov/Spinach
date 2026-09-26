# kernel/optimcon/ctrl_trajan.m

- Signature: `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`

## Purpose

Diagnostic plotting function for optimal control module. Plots trajectory and control pulse analysis. Syntax: ctrl_trajan(spin_system,waveform,trajectory,fidelities)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

## Parameters / inputs

- waveform -waveform, as supplied to a user-end
- function, such as grape_xy.m
- traj_data -a cell array of trajectory data struc-
- tures returned by GRAPE, one per en-
- semble member
- fidelities -fidelities array, as returned by
- user-end functions, such as grape_xy
- Note: this function is called internally by the optimal cont-
- rol module, you should not be calling it directly. All
- settings should be specified in the call to optimcon.m
- when the optimal control problem is set up.

## Implementation structure

- Diagnostic plotting function for optimal control module. Plots
- trajectory and control pulse analysis. Syntax:
- ctrl_trajan(spin_system,waveform,trajectory,fidelities)
- waveform -waveform, as supplied to a user-end
- function, such as grape_xy.m
- traj_data -a cell array of trajectory data struc-
- tures returned by GRAPE, one per en-
- semble member
- fidelities -fidelities array, as returned by
- user-end functions, such as grape_xy
- Note: this function is called internally by the optimal cont-
- rol module, you should not be calling it directly. All
