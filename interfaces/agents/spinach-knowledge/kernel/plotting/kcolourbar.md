# kernel/plotting/kcolourbar.m

- Signature: `kcolourbar(x)`

## Purpose

Creates or updates the colour bar for the current axes and sets its label. Syntax: `kcolourbar(x)`.

## Physical / mathematical content

## Numerical / algorithmic content

- Sets the colour-bar tick-label interpreter to `latex` and its font size to 12; sets the label interpreter to `latex` and its font size to 13.

## Parameters / inputs

- `x` - optional character array used as the label; defaults to `''`. A non-character input raises an error.

## Outputs

- Creates or updates the colour bar in the current axes; returns no output.
