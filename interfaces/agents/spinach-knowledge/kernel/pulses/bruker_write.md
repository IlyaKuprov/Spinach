# kernel/pulses/bruker_write.m

- Signature: `bruker_write(X,Y,dt,file_name)`

## Purpose

Saves pulses in Bruker format. The result is a text file with a list of amplitudes and phases, usable in TopSpin. Syntax: bruker_write(X,Y,dt,file_name)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- X -in-phase channel pulse amplitudes, a
- column vector in Hz
- Y -in-phase channel pulse amplitudes, a
- column vector in Hz
- dt -time slice duration, seconds
- file_name -name of output file with .txt exten-
- sion, a character string

## Outputs

- the function writes an ASCII text file

## Implementation structure

- Saves pulses in Bruker format. The result is a text file with a
- list of amplitudes and phases, usable in TopSpin. Syntax:
- bruker_write(X,Y,dt,file_name)
- X -in-phase channel pulse amplitudes, a
- column vector in Hz
- Y -in-phase channel pulse amplitudes, a
- dt -time slice duration, seconds
- file_name -name of output file with .txt exten-
- sion, a character string
- the function writes an ASCII text file
- Check consistency
- Get amplitudes and phases
