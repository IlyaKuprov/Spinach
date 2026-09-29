# kernel/utilities/min_int_type.m

Source: [kernel/utilities/min_int_type.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/min_int_type.m)

## Purpose

Returns the minimum MATLAB integer data type sufficient to store a specified maximum value. The header comment notes this is useful in many indexing operations in the Spinach kernel where double precision would be a massive overkill.

## Behaviour

- Syntax: `type=min_int_type(max_val,issigned)`.
- The function first runs a consistency check (`grumble`) on the inputs.
- For `'signed'`, it compares `max_val` against `intmax` of `int8`, `int16`, `int32`, and `int64` in ascending order, returning the first type whose maximum covers `max_val`; if `max_val` exceeds `intmax('int64')`, it errors with `Matlab's signed integer types cannot go that far.`
- For `'unsigned'`, it compares `max_val` against `intmax` of `uint8`, `uint16`, `uint32`, and `uint64` in ascending order, returning the first type whose maximum covers `max_val`; if `max_val` exceeds `intmax('uint64')`, it errors with `Matlab's unsigned integer types cannot go that far.`
- Any other `issigned` value triggers the error `unrecognised sign handling type.`
- The consistency check requires `max_val` to be numeric, scalar, real, an integer (`mod(max_val,1)~=0` fails), and at least 1; otherwise it errors with `max_val must be a positive real integer.` It also requires `issigned` to be a character array equal to `'signed'` or `'unsigned'`; otherwise it errors with `valid valued for issigned are 'signed' and 'unsigned'.`

## Inputs and outputs

Inputs:

- `max_val` — maximum value that the integer type must cover; must be a positive real integer scalar.
- `issigned` — whether the integer needs to cover negative values; `'signed'` or `'unsigned'`.

Output:

- `type` — MATLAB data type to use (`'int8'`, `'int16'`, `'int32'`, `'int64'`, `'uint8'`, `'uint16'`, `'uint32'`, or `'uint64'`).

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=min_int_type.m>
