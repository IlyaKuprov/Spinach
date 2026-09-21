# kernel/pulses/restrans.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/restrans.m`
- Signature: `[X,Y,dt]=restrans(X_user,Y_user,dt_user,omega,Q,model,up_factor)`
- Total lines: 196

## Purpose

RLC circuit response calculation -converts a waveform from the ideal shape emitted by the instrument into the shape that comes out of the RLC circuit of the probe. Syntax: [X,Y,dt]=restrans(X_user,Y_user,dt_user,... omega,Q,model,up_factor)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- X_user -in-phase part of the rotating frame
- pulse waveform, a column vector of
- real numbers
- Y_user -out-of-phase part of the rotating
- frame pulse waveform, a column vec-
- tor of real numbers
- dt_user -time slice duration, seconds
- omega -RLC circuit resonance frequency in
- radians per second, a real number
- Q -RLC circuit quality factor, a real
- positive number
- model -input signal model, use 'pwc' for
- piecewise-constant, and 'pwl' for
- piecewise-linear input; time shift
- compensation for piecewise-linear
- is requested by 'pwl_tsc'
- up_factor -the output waveform will have more
- discretisation points than the in-
- put waveform by this factor, about
- 100 is a safe guess

## Outputs

- X -in-phase part of the rotating frame
- pulse waveform distorted by the RLC
- response, a column vector of real
- numbers
- Y -out-of-phase part of the rotating
- frame pulse waveform distorted by
- the RLC response, a column vector
- of real numbers
- dt -slice duration in the distorted wave-
- form, seconds

## Implementation structure

- RLC circuit response calculation -converts a waveform from the
- ideal shape emitted by the instrument into the shape that comes
- out of the RLC circuit of the probe. Syntax:
- [X,Y,dt]=restrans(X_user,Y_user,dt_user,...
- omega,Q,model,up_factor)
- X_user -in-phase part of the rotating frame
- pulse waveform, a column vector of
- real numbers
- Y_user -out-of-phase part of the rotating
- frame pulse waveform, a column vec-
- tor of real numbers
- dt_user -time slice duration, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cumsum()`, `atan2()`, `lsim()`, `lowpass()`, `circuit_time_grid()`, `strcmp()`, `kxlabel()`, `xlim()`, `kylabel()`, `klegend()`, `isscalar()`, `iscolumn()`, `ischar()`, `ismember()`.
