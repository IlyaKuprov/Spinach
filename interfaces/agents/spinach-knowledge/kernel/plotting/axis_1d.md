# kernel/plotting/axis_1d.m

- Signature: `[ax,ax_label]=axis_1d(spin_system,parameters)`

## Purpose

Generates axis ticks for plotting 1D spectra. Syntax: [ax,ax_label]=axis_1d(spin_system,parameters)

## Physical / mathematical content

A one-dimensional spectrum axis expressed in the requested frequency, field, g-tensor, or digitisation units.

## Numerical / algorithmic content

The function constructs the Hz axis with `ft_axis(offset,sweep,zerofill)` for a sweep width (requiring `offset`), or with `linspace(sweep(1),sweep(2),zerofill)` for two increasing endpoints (with no `offset`). It converts the resulting axis to the requested units and returns the axis and plot label; it does not perform an FFT or apodisation.

## Parameters / inputs

- `sweep` — one real numeric value for sweep width in Hz, or two increasing real numeric endpoints in Hz.
- `zerofill` — positive integer number of points in the axis.
- `offset` — spectrum centre offset relative to the magnet frequency, in Hz; required for a one-value sweep and not supplied for a two-endpoint sweep.
- `axis_units` — one of `ppm`, `Gauss`, `mT`, `Hz`, `kHz`, `MHz`, `MHz-labframe`, `GHz`, `GHz-labframe`, `gtensor`, or `points`.
- `spins` — one-element cell array naming an isotope present in the spin system, such as `{'1H'}`.

## Outputs

- `ax` — row vector of axis values.
- `ax_label` — label for the selected axis units and spin.

Magnetic-field units use the free-electron g-tensor for conversion.
