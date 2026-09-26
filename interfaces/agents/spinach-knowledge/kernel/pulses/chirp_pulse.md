# kernel/pulses/chirp_pulse.m

- Signature: `[Cx,Cy,durs,ints,amps,phis,frqs]=...`

## Purpose

Chirp pulse waveform with a sine bell power or a quarter-sine amplitude fade-in and fade-out. Generates unidirectional chir- ps or saltire chirps which are super-positions of two counter- sweeping chirps. Syntax: [Cx,Cy,durs,ints,amps,phis,frqs]=... chirp_pulse(npts,dur,bwidth,smp,type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- npts -number of discretization points in
- the waveform
- duration -pulse duration, seconds
- bwidth -chirp sweep bandwidth around
- zero frequency, Hz
- type -'wurst', 'smoothed', or 'saltire'; the
- default is uniform time grid, to get
- adaptive sampling, add '-adaptive'
- smp -smoothing parameter; for 'wurst', this
- is the power in
- 1-|sin(x)^smp|
- as x approaches pi/2 at either the edge
- of the pulse. For 'smoothed' and 'salti-
- re', this is the fraction of the pulse
- duration (in percent) that is affected
- by a sine bell fade-in and fade-out: 0
- means square amplitude envelope and 50
- means sine bell envelope.

## Outputs

- Cx -real part of the waveform, calibrated to
- produce an inversion pulse, rad/s
- Cy -imag part of the waveform, calibrated to
- produce an inversion pulse, rad/s
- durs -slice durations for piecewise-constant
- approximation, seconds
- ints -interval durations for piecewise-linear
- approximation, seconds
- amps -waveform amplitudes, rad/s
- phis -waveform phases, rad
- frqs -waveform frequencies, Hz
- intv_grid -normalised interval grid, npts-1 elements
- Note: Cy is zero for the saltire pulse, this radically changes
- its phase and amplitude profiles.

## Implementation structure

- Chirp pulse waveform with a sine bell power or a quarter-sine
- amplitude fade-in and fade-out. Generates unidirectional chir-
- ps or saltire chirps which are super-positions of two counter-
- sweeping chirps. Syntax:
- [Cx,Cy,durs,ints,amps,phis,frqs]=...
- chirp_pulse(npts,dur,bwidth,smp,type)
- npts -number of discretization points in
- the waveform
- duration -pulse duration, seconds
- bwidth -chirp sweep bandwidth around
- zero frequency, Hz
- type -'wurst', 'smoothed', or 'saltire'; the
