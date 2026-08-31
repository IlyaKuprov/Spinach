# kernel/pulses/bloch_siegert.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/bloch_siegert.m`
- Signature: `[ctrl_opers,...`
- Total lines: 114

## Purpose

Applies Bloch-Siegert corrections to Cartesian control pulses. Takes control operators and their amplitude arrays and returns augmented arrays with the Bloch-Siegert response operator of each control chan- nel appended as a virtual channel whose coefficient vector is the square of the corresponding control amplitude vector. The output arguments are ready for the shaped_pulse_xy function, which then propagates the sys

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Check consistency; implemented by `grumble(spin_system,ctrl_opers,ctrl_coefs)`.
- Lines 52-53: Count the physical channels; implemented by `n_ctrls=numel(ctrl_opers)`.
- Lines 55-56: Append the response operator channels; implemented by `for n=1:n_ctrls`.

### Control flow inferred from the code

- Line 56: `for` loop over `n=1:n_ctrls`.

### Key state/data transformations

- Lines 46-47: computes `ctrl_coefs]` using `ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers, ctrl_coefs)`.
- Lines 53: computes `n_ctrls` using `n_ctrls=numel(ctrl_opers)`.
- Lines 57: computes `ctrl_opers{n_ctrls+n}` using `ctrl_opers{n_ctrls+n}=spin_system.control.resp_ops{n}`.
- Lines 58: computes `ctrl_coefs{n_ctrls+n}` using `ctrl_coefs{n_ctrls+n}=ctrl_coefs{n}.^2`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(spin_system,ctrl_opers,ctrl_coefs)`.
  - Representative operation: `if ~isstruct(spin_system)`.
  - Representative operation: `error('spin_system must be a structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system object containing the
- Bloch-Siegert settings added by optimcon()
- ctrl_opers -a cell array of control operators, one per
- control channel
- ctrl_coefs -a cell array of control coefficient vectors
- in rad/s, one vector per control channel

## Outputs

- ctrl_opers -an augmented cell array of control opera-
- tors, with the response operator of each
- channel appended
- ctrl_coefs -an augmented cell array of coefficient vec-
- tors, with the squared amplitude vector of
- each channel appended
- Note: the augmented arrays are only valid for the piecewise-constant
- propagation methods of shaped_pulse_xy; piecewise-linear methods
- would interpolate the squared coefficients, which does not cor-
- respond to the square of the interpolated control amplitude.

## Implementation structure

- Applies Bloch-Siegert corrections to Cartesian control pulses. Takes
- control operators and their amplitude arrays and returns augmented
- arrays with the Bloch-Siegert response operator of each control chan-
- nel appended as a virtual channel whose coefficient vector is the
- square of the corresponding control amplitude vector. The output
- arguments are ready for the shaped_pulse_xy function, which then
- propagates the system under the same slice generators that the GRAPE
- engines use when Bloch-Siegert corrections are enabled. Syntax:
- [ctrl_opers,...
- ctrl_coefs]=bloch_siegert(spin_system,ctrl_opers,...
- ctrl_coefs)
- spin_system -Spinach spin system object containing the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bloch_siegert()`, `grumble()`, `isstruct()`, `isfield()`, `optimcon()`, `islogical()`, `isscalar()`, `iscell()`, `isvector()`.
