# kernel/pulses/vg_pulse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/vg_pulse.m`
- Signature: `waveform=vg_pulse(pulse_name,npoints,duration)`
- Total lines: 136

## Purpose

Veshtort-Griffin shaped pulses, generated from tables given in There are good reasons to believe (see Section 2.2 of the paper) that these are the best possible pulses within their design specifications and basis sets. Syntax: waveform=vg_pulse(pulse_name,npoints,duration)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(pulse_name,npoints,duration)`.
- Lines 37-39: Stockpile names; implemented by `names={'E0A' 'E0B' 'E100A' 'E100B' 'E200A' 'E200D' 'E200F' 'E300C' 'E300F' 'E400B' 'E300A' 'E500A' 'E500B' 'E500C' 'E600A' 'E600C' 'E600F' 'E800A' 'E800B' 'E1000B'}`.
- Lines 41-43: Stockpile cosine coefficients; implemented by `cosines=[ -0.763 -0.763 -0.734 0.265 0.282 0.276 0.251 0.228 0.269 0.231 0.222 0.240 0.242 0.222 -0.739 0.239 0.249 0.238 0.235 0.254`.
- Lines 65-67: Stockpile sine coefficients; implemented by `sines=[ 0.000 0.000 0.092 -2.104 1.459 -1.270 -1.273 -1.279 -0.788 -0.533 -0.772 -0.681 -3.105 -0.215 -0.360 0.249 -2.489 -0.238 -0.773 0.242`.
- Lines 88-89: Parse the pulse name; implemented by `pulse_number=find(strcmp(names,pulse_name))`.
- Lines 92-93: Preallocate the pulse; implemented by `waveform=zeros(1,npoints)`.
- Lines 95-96: Build the time grid; implemented by `time_grid=linspace(0,2*pi,npoints)`.
- Lines 98-99: Compute cosine terms; implemented by `for k=0:20`.
- Lines 103-104: Compute sine terms; implemented by `for k=1:20`.
- Lines 108-109: Normalise the waveform; implemented by `waveform=2*pi*waveform/duration`.

### Control flow inferred from the code

- Line 90: conditional branch on `isempty(pulse_number), error('incorrect pulse name.'); end`.
- Line 99: `for` loop over `k=0:20`.
- Line 104: `for` loop over `k=1:20`.

### Key state/data transformations

- Lines 38-39: computes `names` using `names={'E0A' 'E0B' 'E100A' 'E100B' 'E200A' 'E200D' 'E200F' 'E300C' 'E300F' 'E400B' 'E300A' 'E500A' 'E500B' 'E500C' 'E600A' 'E600C' 'E600F' 'E800A' 'E800B' 'E1000B'}`.
- Lines 42-43: computes `cosines` using `cosines=[ -0.763 -0.763 -0.734 0.265 0.282 0.276 0.251 0.228 0.269 0.231 0.222 0.240 0.242 0.222 -0.739 0.239 0.249 0.238 0.235 0.254`.
- Lines 66-67: computes `sines` using `sines=[ 0.000 0.000 0.092 -2.104 1.459 -1.270 -1.273 -1.279 -0.788 -0.533 -0.772 -0.681 -3.105 -0.215 -0.360 0.249 -2.489 -0.238 -0.773 0.242`.
- Lines 89: computes `pulse_number` using `pulse_number=find(strcmp(names,pulse_name))`.
- Lines 93: computes `waveform` using `waveform=zeros(1,npoints)`.
- Lines 96: computes `time_grid` using `time_grid=linspace(0,2*pi,npoints)`.

### Local helper functions

- Line 114: `grumble()` — `function grumble(pulse_name,npoints,duration)`.
  - Representative operation: `if ~ischar(pulse_name)`.
  - Representative operation: `error('pulse_name must be a character string.')`.

## Parameters / inputs

- pulse_name -a character string, one of the following: E0A,
- E0B, E100A, E100B, E200A, E200D, E200F, E300C,
- E300F, E400B, E300A, E500A, E500B, E500C, E600A,
- E600C, E600F, E800A, E800B, E1000B
- npoints -number of discrete time intervals in the pulse
- duration -duration of the pulse, seconds

## Outputs

- waveform -amplitude of the pulse at each interval (there
- is no phase modulation), normalised to produce
- a 90-degree pulse, rad/s

## Implementation structure

- Veshtort-Griffin shaped pulses, generated from tables given in
- There are good reasons to believe (see Section 2.2 of the paper) that
- these are the best possible pulses within their design specifications
- and basis sets. Syntax:
- waveform=vg_pulse(pulse_name,npoints,duration)
- pulse_name -a character string, one of the following: E0A,
- E0B, E100A, E100B, E200A, E200D, E200F, E300C,
- E300F, E400B, E300A, E500A, E500B, E500C, E600A,
- E600C, E600F, E800A, E800B, E1000B
- npoints -number of discrete time intervals in the pulse
- duration -duration of the pulse, seconds
- waveform -amplitude of the pulse at each interval (there

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmp()`, `cosines()`, `sines()`, `ischar()`.
