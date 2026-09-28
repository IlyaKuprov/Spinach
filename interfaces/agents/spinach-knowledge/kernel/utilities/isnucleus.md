# kernel/utilities/isnucleus.m

- Signature: `verdict=isnucleus(spin_spec)`

## Purpose

Classifies a valid Spinach particle specification as a nucleus or not. The function first requires a character string and validates the specification with `spin`.

## Parameters / inputs

- `spin_spec` - character string containing a Spinach particle specification.

## Outputs

- `verdict` - true when the specification passes the function's nucleus test; false otherwise.

## Numerical / algorithmic content

The function returns false when the first character is one of `E`, `C`, `V`, or `T`, or when the complete specification is one of `G`, `E`, `N`, or `M`. It returns true otherwise.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=isnucleus.m)
