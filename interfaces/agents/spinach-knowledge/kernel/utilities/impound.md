# kernel/utilities/impound.m

- Signature: `answer=impound(varargin)`

## Purpose

Returns all arguments received by the function, collected in a cell array. This is useful for passing information back from Spinach wrappers when `impound` is used as a pulse sequence.

## Parameters / inputs

- `varargin` - any number of input parameters of any type.

## Outputs

- `answer` - a cell array containing the input parameters, in their original order.

## Implementation

The function assigns `varargin` directly to `answer`; it does not transform or interpret the inputs.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=impound.m)
