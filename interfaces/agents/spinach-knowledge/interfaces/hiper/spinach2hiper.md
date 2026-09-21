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
