# kernel/utilities/iselectron.m

- Signature: `verdict=iselectron(spin_spec)`

## Purpose

Tests whether a Spinach particle specification denotes an electron. The function checks that `spin_spec` is a character string and a valid Spinach particle specification, then returns true when its first character is `E` and false otherwise.

## Parameters / inputs

- `spin_spec` - a character string containing a Spinach particle specification.

## Outputs

- `verdict` - true if the first character of `spin_spec` is `E`, false otherwise.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=iselectron.m)
