# kernel/pulses/bruker_write.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/bruker_write.m`
- Signature: `bruker_write(X,Y,dt,file_name)`
- Total lines: 108

## Purpose

Saves pulses in Bruker format. The result is a text file with a list of amplitudes and phases, usable in TopSpin. Syntax: bruker_write(X,Y,dt,file_name)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cartesian2polar()`, `rad2deg()`, `wrapTo2Pi()`, `datetime()`, `string()`, `num2str()`, `writelines()`, `writematrix()`, `iscolumn()`, `isscalar()`, `ischar()`.
