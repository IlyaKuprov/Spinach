# kernel/utilities/impound.m

## Purpose

`impound.m` packages everything it receives into a cell array and returns it back. It is useful for pulling information back from various Spinach wrappers by calling it as a pulse sequence.

## Behaviour

The function is defined as `answer=impound(varargin)`. It returns what was received by assigning `answer=varargin`, so all input arguments are collected into a single cell array. The function contains no other logic.

## Inputs and outputs

Inputs:

- `varargin` — any number of parameters of any type.

Outputs:

- `answer` — all input parameters as a cell array.

## References

- Source: [kernel/utilities/impound.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/impound.m)
- Spinach Wiki: [impound.m](https://spindynamics.org/wiki/index.php?title=impound.m)
