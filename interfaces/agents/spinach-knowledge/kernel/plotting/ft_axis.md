# kernel/plotting/ft_axis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/ft_axis.m`
- Signature: `ax=ft_axis(offset,sweep,npoints)`
- Total lines: 65

## Purpose

Fourier transform axis ticks generator that accounts for the periodicity and correctly folds the edge frequency. Syntax: ax=ft_axis(offset,sweep,npoints)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- offset -centre frequency
- sweep -frequency range
- npoints -number of points

## Outputs

- ax -row vector of axis ticks

## Implementation structure

- Fourier transform axis ticks generator that accounts for the
- periodicity and correctly folds the edge frequency. Syntax:
- ax=ft_axis(offset,sweep,npoints)
- offset -centre frequency
- sweep -frequency range
- npoints -number of points
- ax -row vector of axis ticks
- Check consistency
- Axis with an extra point
- Odd and even point counts
- If odd, drop and shift
- If even, just drop

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
