# kernel/plotting/sweep2ticks.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/sweep2ticks.m`
- Signature: `axis_hz=sweep2ticks(offs,sweep,npoints)`
- Total lines: 55

## Purpose

Converts offset-sweep-npoints specification into axis ticks in Hz. The function returns the frequency axis of the spectrum, suitable for use in Matlab functions like plot(). Syntax: axis_hz=sweep2ticks(offs,sweep,npoints)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- offs -offset from carrier frequency, Hz
- sweep -sweep width, Hz
- npoints -number of points in the spectrum

## Outputs

- axis_hz -a column vector of axis ticks, Hz

## Implementation structure

- Converts offset-sweep-npoints specification into axis ticks in Hz.
- The function returns the frequency axis of the spectrum, suitable
- for use in Matlab functions like plot(). Syntax:
- axis_hz=sweep2ticks(offs,sweep,npoints)
- offs -offset from carrier frequency, Hz
- sweep -sweep width, Hz
- npoints -number of points in the spectrum
- axis_hz -a column vector of axis ticks, Hz
- Check consistency
- Build the axis
- Consistency enforcement
- Spinach code is clear, useful and elegant because the program is the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
