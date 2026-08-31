# interfaces/hiper/spinach2hiper.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/hiper/spinach2hiper.m`
- Signature: `spinach2hiper(file_name,amp,phi,off,dt)`
- Total lines: 83

## Purpose

Exports phase-modulated optimal control waveforms into the format expected by Graham Smith's HiPER instrument. Syntax: spinach2hiper(file_name,amp,phi,off,dt)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(file_name,amp,phi,off,dt)`.
- Lines 36-37: Convert and wrap phases; implemented by `phi=wrapTo360(180*phi(:)/pi)`.
- Lines 39-40: Slice timing and offsets; implemented by `off=1e-6*repmat(off,size(phi))`.
- Lines 44-46: Make a table for export; implemented by `pulse_table=table(times,off,phi,amp(:),'VariableNames', {'time_ns','freq_MHz','phase_deg','amplitude'})`.
- Lines 48-49: Export the table into a CSV file; implemented by `writetable(pulse_table,[file_name '.csv'])`.

### Key state/data transformations

- Lines 37: computes `phi` using `phi=wrapTo360(180*phi(:)/pi)`.
- Lines 40: computes `off` using `off=1e-6*repmat(off,size(phi))`.
- Lines 41: computes `dt` using `dt=1e9*repmat(dt,size(phi))`.
- Lines 42: computes `times` using `times=cumsum(dt)-dt`.
- Lines 45-46: computes `pulse_table` using `pulse_table=table(times,off,phi,amp(:),'VariableNames', {'time_ns','freq_MHz','phase_deg','amplitude'})`.

### Local helper functions

- Line 54: `grumble()` — `function grumble(file_name,amp,phi,off,dt)`.
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

## Parameters / inputs

- file_name -CSV file name, a character string
- without the extension
- amp -a vector of amplitudes
- phi -a vector of phases in radians
- off -transmitter offset in Hz
- dt -waveform slice duration, seconds

## Outputs

- this function writes a file
- Note: in practice, try both positive and negative pha-
- ses -some instruments count phases clockwise,
- others counterclockwise.

## Implementation structure

- Exports phase-modulated optimal control waveforms into the
- format expected by Graham Smith's HiPER instrument. Syntax:
- spinach2hiper(file_name,amp,phi,off,dt)
- file_name -CSV file name, a character string
- without the extension
- amp -a vector of amplitudes
- phi -a vector of phases in radians
- off -transmitter offset in Hz
- dt -waveform slice duration, seconds
- this function writes a file
- Note: in practice, try both positive and negative pha-
- ses -some instruments count phases clockwise,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `wrapTo360()`, `phi()`, `cumsum()`, `table()`, `amp()`, `writetable()`, `ischar()`, `isvector()`, `any()`, `isscalar()`.
