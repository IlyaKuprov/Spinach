# kernel/optimcon/ctrl_trajan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/ctrl_trajan.m`
- Signature: `ctrl_trajan(spin_system,waveform,traj_data,fidelities)`
- Total lines: 791

## Purpose

Diagnostic plotting function for optimal control module. Plots trajectory and control pulse analysis. Syntax: ctrl_trajan(spin_system,waveform,trajectory,fidelities)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `setdiff()`, `ismember()`, `set()`, `scale_figure()`, `subplot()`, `waveform()`, `fliplr()`, `spectrogram()`, `atan2()`, `cat()`, `image()`, `hsv2rgb()`, `ktitle()`, `num2str()`, `kylabel()`.
