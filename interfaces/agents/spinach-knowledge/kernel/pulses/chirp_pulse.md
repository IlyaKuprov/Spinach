# kernel/pulses/chirp_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/chirp_pulse.m`
- Signature: `[Cx,Cy,durs,ints,amps,phis,frqs]=...`
- Total lines: 203

## Purpose

Chirp pulse waveform with a sine bell power or a quarter-sine amplitude fade-in and fade-out. Generates unidirectional chir- ps or saltire chirps which are super-positions of two counter- sweeping chirps. Syntax: [Cx,Cy,durs,ints,amps,phis,frqs]=... chirp_pulse(npts,dur,bwidth,smp,type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 72-73: Check consistency; implemented by `grumble(npts,dur,bwidth,smp,type)`.
- Lines 75-76: Recursive call for saltire; implemented by `if contains(type,'saltire')`.
- Lines 78-79: Saltire pulses are not monochromatic; implemented by `if nargout>6, error('saltire pulse is not monochromatic.'); end`.
- Lines 81-82: Request a smoothed chirp; implemented by `type=replace(type,'saltire','smoothed')`.
- Lines 85-86: Zero out Y component and update parameters; implemented by `Cy=zeros(size(Cy)); amps=abs(Cx); phis=pi*(Cx<0); return`.
- Lines 90-91: Decide the time grid; implemented by `if contains(type,'-adaptive')`.
- Lines 93-94: Adaptive normalised time grid; implemented by `time_grid=linspace(-1.0,1.0,npts)`.
- Lines 98-99: For PWL: interval durations; implemented by `ints=dur*diff(time_grid)`.
- Lines 101-102: For PWC: slice durations; implemented by `aux_grid=linspace(-1.0,1.0,npts+1)`.
- Lines 108-109: Uniform normalised time grid; implemented by `time_grid=linspace(-0.5,0.5,npts)`.
- Lines 114-115: For PWC: slice durations; implemented by `durs=(dur/npts)*ones(1,npts)`.
- Lines 119-120: Phase and frequency sequence; implemented by `phis=pi*dur*bwidth*(time_grid.^2)`.
- Lines 123-124: Sampling adequacy check; implemented by `freq_grid=linspace(-bwidth/2,bwidth/2,npts-1)`.
- Lines 130-131: Amplitude sequence; implemented by `switch type`.
- Lines 135-136: Pulse amplitude envelope; implemented by `amps=1-abs(sin(pi*time_grid).^smp)`.
- Lines 140-141: Number of points to fade either side; implemented by `np_do_smooth=nnz(time_grid>(50-smp)/100)`.
- Lines 143-144: Number of points to leave intact; implemented by `np_no_smooth=npts-2*np_do_smooth`.
- Lines 146-147: Compute the amplitude envelope; implemented by `fade_in=sin(linspace(0,pi/2,np_do_smooth))`.

### Control flow inferred from the code

- Line 76: conditional branch on `contains(type,'saltire')`.
- Line 79: conditional branch on `nargout>6, error('saltire pulse is not monochromatic.'); end`.
- Line 91: conditional branch on `contains(type,'-adaptive')`.
- Line 126: conditional branch on `any(abs(phi_jumps)>pi,'all')&&(nargout<7)`.
- Line 131: dispatches on `type`; cases `{'wurst','wurst-adaptive'}`, `{'smoothed','smoothed-adaptive'}`.

### Key state/data transformations

- Lines 82: computes `type` using `type=replace(type,'saltire','smoothed')`.
- Lines 83: computes `[Cx,Cy,durs,ints]` using `[Cx,Cy,durs,ints]=chirp_pulse(npts,dur,bwidth,smp,type)`.
- Lines 86: computes `Cy` using `Cy=zeros(size(Cy)); amps=abs(Cx); phis=pi*(Cx<0); return`.
- Lines 94: computes `time_grid` using `time_grid=linspace(-1.0,1.0,npts)`.
- Lines 99: computes `ints` using `ints=dur*diff(time_grid)`.
- Lines 102: computes `aux_grid` using `aux_grid=linspace(-1.0,1.0,npts+1)`.
- Lines 115: computes `durs` using `durs=(dur/npts)*ones(1,npts)`.
- Lines 120: computes `phis` using `phis=pi*dur*bwidth*(time_grid.^2)`.
- Lines 121: computes `frqs` using `frqs=bwidth*time_grid`.
- Lines 124: computes `freq_grid` using `freq_grid=linspace(-bwidth/2,bwidth/2,npts-1)`.
- Lines 125: computes `phi_jumps` using `phi_jumps=2*pi*abs(freq_grid).*ints`.
- Lines 136: computes `amps` using `amps=1-abs(sin(pi*time_grid).^smp)`.
- Lines 141: computes `np_do_smooth` using `np_do_smooth=nnz(time_grid>(50-smp)/100)`.
- Lines 144: computes `np_no_smooth` using `np_no_smooth=npts-2*np_do_smooth`.
- Lines 147: computes `fade_in` using `fade_in=sin(linspace(0,pi/2,np_do_smooth))`.
- Lines 148: computes `fade_out` using `fade_out=fliplr(fade_in)`.
- Lines 162: computes `[Cx,Cy]` using `[Cx,Cy]=polar2cartesian(amps,phis)`.

### Local helper functions

- Line 167: `grumble()` — `function grumble(npts,dur,bwidth,smp,type)`.
  - Representative operation: `if (~isnumeric(dur))||(~isreal(dur))|| (numel(dur)~=1)||(~isfinite(dur))||(dur<=0)`.
  - Representative operation: `(numel(dur)~=1)||(~isfinite(dur))||(dur<=0)`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `chirp_pulse()`, `grumble()`, `contains()`, `replace()`, `sign()`, `diff()`, `any()`, `nnz()`, `fliplr()`, `polar2cartesian()`, `ischar()`, `ismember()`.
