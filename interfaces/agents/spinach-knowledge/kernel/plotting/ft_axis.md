# kernel/plotting/ft_axis.m

- Signature: `ax=ft_axis(offset,sweep,npoints)`

## Purpose

Fourier transform axis ticks generator that accounts for the periodicity and correctly folds the edge frequency. Syntax: ax=ft_axis(offset,sweep,npoints)

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

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
