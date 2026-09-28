# kernel/overloads/@struct/mtimes.m

- Signature: `str_out=mtimes(M,str_in)`

## Purpose

Multiply every numeric leaf field of a structure by `M`, recursively processing nested structures.

## Parameters / inputs

- `M` - numeric multiplier.
- `str_in` - structure whose fields are numeric values or nested structures of such values.

## Outputs

- `str_out` - structure with the same field layout and each numeric leaf replaced by `M*value`.

## Implementation structure

The function checks that `M` is numeric and `str_in` is a structure, then visits each field and applies `M*field`. MATLAB's matrix multiplication rules apply to each leaf.
