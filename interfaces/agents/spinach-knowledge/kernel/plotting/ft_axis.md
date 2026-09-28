# kernel/plotting/ft_axis.m

- Signature: `ax=ft_axis(offset,sweep,npoints)`

## Purpose

Generates a row vector of frequency-axis ticks centred on `offset` over the range `sweep`, with spacing `sweep/npoints`.

## Physical / mathematical content

- The ticks span a frequency interval of width `sweep` around `offset`. For odd `npoints`, the first endpoint is dropped and the remaining ticks are shifted by half a spacing; for even `npoints`, the last endpoint is dropped.

## Numerical / algorithmic content

- The function starts with `npoints+1` equally spaced values over `[-sweep/2,sweep/2]`, applies the parity-dependent endpoint adjustment, then adds `offset`.
- `offset` must be a real numeric scalar, `sweep` a positive real numeric scalar, and `npoints` a real numeric integer greater than 2.

## Parameters / inputs

- offset -centre frequency
- sweep -frequency range
- npoints -number of points

## Outputs

- ax -row vector of axis ticks
