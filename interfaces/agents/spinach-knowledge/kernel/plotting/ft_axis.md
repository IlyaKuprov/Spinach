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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(offset,sweep,npoints)`.
- Lines 27-28: Axis with an extra point; implemented by `ax=linspace(-sweep/2,sweep/2,npoints+1)`.
- Lines 30-31: Odd and even point counts; implemented by `if mod(npoints,2)==1`.
- Lines 33-34: If odd, drop and shift; implemented by `ax=ax(2:end)-(ax(2)-ax(1))/2+offset`.
- Lines 38-39: If even, just drop; implemented by `ax=ax(1:(end-1))+offset`.

### Control flow inferred from the code

- Line 31: conditional branch on `mod(npoints,2)==1`.

### Key state/data transformations

- Lines 28: computes `ax` using `ax=linspace(-sweep/2,sweep/2,npoints+1)`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(offset,sweep,npoints)`.
  - Representative operation: `if (~isnumeric(offset))||(~isreal(offset))|| (~isscalar(offset))`.
  - Representative operation: `(~isscalar(offset))`.

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
