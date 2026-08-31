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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 59-60: Check consistency; implemented by `grumble(X_user,Y_user,dt_user,omega,Q,model,up_factor)`.
- Lines 62-63: 16 x Nyquist oversampling; implemented by `dt=pi/(16*omega)`.
- Lines 65-66: Signal model; implemented by `switch model`.
- Lines 68-69: Piecewise-constant; implemented by `case 'pwc'`.
- Lines 71-72: Assume that the user has slice midpoints; implemented by `nslices=numel(X_user); tmax=nslices*dt_user`.
- Lines 75-76: Time grid required by the RLC circuit; implemented by `circuit_time_grid=linspace(0,tmax,tmax/dt+1)'`.
- Lines 78-79: Interpolate user waveform onto RLC circuit time grid; implemented by `X0=interp1(user_time_grid,X_user,circuit_time_grid,'nearest','extrap')`.
- Lines 82-83: Piecewise-linear; implemented by `case {'pwl','pwl_tsc'}`.
- Lines 85-86: Assume that the user has slice edges; implemented by `nslices=numel(X_user)-1; tmax=nslices*dt_user`.
- Lines 92-93: Interpolate user waveform onto RLC circuit time grid; implemented by `X0=interp1(user_time_grid,X_user,circuit_time_grid,'linear')`.
- Lines 98-99: Complain and bomb out; implemented by `error('unknown signal model.')`.
- Lines 103-104: Fancy colour palette; implemented by `rgb_rd=[0.8500 0.3250 0.0980]`.
- Lines 109-110: Diagnostic plotting if no outputs; implemented by `if nargout==0`.
- Lines 117-118: Generate wall clock input signal; implemented by `inp_amp=sqrt(X0.^2+Y0.^2); inp_phi=atan2(Y0,X0)`.
- Lines 121-122: Build the RLC circuit response kernel; implemented by `sys=tf(1/Q,[1/(omega^2) 1/(omega*Q) 1])`.
- Lines 124-125: Apply the RLC circuit response kernel; implemented by `out_signal=lsim(sys,inp_signal,circuit_time_grid)`.
- Lines 127-129: Heterodyne out the carrier frequency; implemented by `X=2*lowpass(out_signal.*sin(omega*circuit_time_grid), 1,64,ImpulseResponse="iir",Steepness=0.95)`.
- Lines 133-134: Downsample the rotating frame waveform; implemented by `down_factor=numel(X)/(up_factor*numel(X_user))`.

### Control flow inferred from the code

- Line 66: dispatches on `model`; cases `'pwc'`, `{'pwl','pwl_tsc'}`.
- Line 110: conditional branch on `nargout==0`.
- Line 140: conditional branch on `strcmp(model,'pwl_tsc')`.
- Line 145: conditional branch on `nargout==0`.

### Key state/data transformations

- Lines 63: computes `dt` using `dt=pi/(16*omega)`.
- Lines 72: computes `nslices` using `nslices=numel(X_user); tmax=nslices*dt_user`.
- Lines 73: computes `user_time_grid` using `user_time_grid=dt_user*(cumsum(ones(nslices,1))-1/2)`.
- Lines 76: computes `circuit_time_grid` using `circuit_time_grid=linspace(0,tmax,tmax/dt+1)'`.
- Lines 79: computes `X0` using `X0=interp1(user_time_grid,X_user,circuit_time_grid,'nearest','extrap')`.
- Lines 80: computes `Y0` using `Y0=interp1(user_time_grid,Y_user,circuit_time_grid,'nearest','extrap')`.
- Lines 104: computes `rgb_rd` using `rgb_rd=[0.8500 0.3250 0.0980]`.
- Lines 105: computes `rgb_bd` using `rgb_bd=[0.0000 0.4470 0.7410]`.
- Lines 106: computes `rgb_rb` using `rgb_rb=rgb_rd+1; rgb_rb=rgb_rb/max(rgb_rb)`.
- Lines 107: computes `rgb_bb` using `rgb_bb=rgb_bd+1; rgb_bb=rgb_bb/max(rgb_bb)`.
- Lines 118: computes `inp_amp` using `inp_amp=sqrt(X0.^2+Y0.^2); inp_phi=atan2(Y0,X0)`.
- Lines 119: computes `inp_signal` using `inp_signal=inp_amp.*cos(omega*circuit_time_grid+inp_phi)`.
- Lines 122: computes `sys` using `sys=tf(1/Q,[1/(omega^2) 1/(omega*Q) 1])`.
- Lines 125: computes `out_signal` using `out_signal=lsim(sys,inp_signal,circuit_time_grid)`.
- Lines 128-129: computes `X` using `X=2*lowpass(out_signal.*sin(omega*circuit_time_grid), 1,64,ImpulseResponse="iir",Steepness=0.95)`.
- Lines 129: computes `1,64,ImpulseResponse` using `1,64,ImpulseResponse="iir",Steepness=0.95)`.
- Lines 130-131: computes `Y` using `Y=2*lowpass(out_signal.*cos(omega*circuit_time_grid), 1,64,ImpulseResponse="iir",Steepness=0.95)`.
- Lines 134: computes `down_factor` using `down_factor=numel(X)/(up_factor*numel(X_user))`.

### Local helper functions

- Line 159: `grumble()` — `function grumble(X_user,Y_user,dt_user,omega,Q,model,up_factor)`.
  - Representative operation: `if (~isnumeric(dt_user))||(~isreal(dt_user))|| (~isscalar(dt_user))||(~isfinite(dt_user))||(dt_user<=0)`.
  - Representative operation: `(~isscalar(dt_user))||(~isfinite(dt_user))||(dt_user<=0)`.

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
