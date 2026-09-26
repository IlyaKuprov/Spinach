# interfaces/hiper/spinach2hiper.m

- Signature: `spinach2hiper(file_name,amp,phi,off,dt)`

## Purpose

Exports phase-modulated optimal control waveforms into the format expected by Graham Smith's HiPER instrument. Syntax: spinach2hiper(file_name,amp,phi,off,dt)

## Physical / mathematical content

## Numerical / algorithmic content

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
