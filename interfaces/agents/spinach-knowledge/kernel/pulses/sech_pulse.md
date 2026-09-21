# kernel/pulses/sech_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/sech_pulse.m`
- Signature: `[Cx,Cy,time_grid,amps,phis]=...`
- Total lines: 93

## Purpose

Hyperbolic secant pulse in Cartesian and amplitude-phase representation. Syntax: [Cx,Cy,time_grid,amps,phis]=... sech_pulse(peak_ampl,freq_mod,phase_mod,dur,npts)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Parameters / inputs

- peak_amp -peak amplitude, rad/s
- freq_mod -frequency modulation parameter, rad/s
- phase_mod -phase modulation parameter, dimless
- dur -pulse duration, seconds
- npts -number of digitisation points

## Outputs

- Cx -a vector of coefficients in front of Sx
- spin operator at each time slice, rad/s
- Cy -a vector of coefficients in front of Sy
- spin operator at each time slice, rad/s
- time_grid -a vector of time grid points, seconds
- amps -a vector of pulse amplitudes at each ti-
- me slice, rad/s
- phis -a vector of pulse phases at each time
- slice (phi=0 at t=0 in the centre), rad
- Example:
- [Cx,Cy,time_grid]=sech_pulse(1,672,5,10.24e-3,1000);
- plot(time_grid,[Cx; Cy]); kgrid; xlim tight;
- kxlabel('time, seconds'); kylabel('amplitude, rad/s');

## Implementation structure

- Hyperbolic secant pulse in Cartesian and amplitude-phase
- representation. Syntax:
- [Cx,Cy,time_grid,amps,phis]=...
- sech_pulse(peak_ampl,freq_mod,phase_mod,dur,npts)
- peak_amp -peak amplitude, rad/s
- freq_mod -frequency modulation parameter, rad/s
- phase_mod -phase modulation parameter, dimless
- dur -pulse duration, seconds
- npts -number of digitisation points
- Cx -a vector of coefficients in front of Sx
- spin operator at each time slice, rad/s
- Cy -a vector of coefficients in front of Sy

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sech_pulse()`, `grumble()`, `sech()`, `cosh()`, `polar2cartesian()`, `isscalar()`.
