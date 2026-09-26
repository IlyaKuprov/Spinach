# kernel/plotting/sweep2ticks.m

- Signature: `axis_hz=sweep2ticks(offs,sweep,npoints)`

## Purpose

Converts offset-sweep-npoints specification into axis ticks in Hz. The function returns the frequency axis of the spectrum, suitable for use in Matlab functions like plot(). Syntax: axis_hz=sweep2ticks(offs,sweep,npoints)

## Physical / mathematical content

## Numerical / algorithmic content

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
