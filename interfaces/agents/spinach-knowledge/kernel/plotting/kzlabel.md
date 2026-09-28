# kernel/plotting/kzlabel.m

- Signature: `kzlabel(varargin)`

## Purpose

House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kzlabel(varargin)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `varargin` — same arguments as those accepted by Matlab's `zlabel` function.

## Outputs

- creates or updates the current axis system

## Implementation structure

- Calls `zlabel` with the supplied arguments and the `latex` interpreter.
- Sets the current axes tick-label interpreter to `latex` and its font size to 12.

[Source reference](https://spindynamics.org/wiki/index.php?title=kzlabel.m)
