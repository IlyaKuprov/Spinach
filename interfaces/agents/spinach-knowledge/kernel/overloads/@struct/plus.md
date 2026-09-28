# kernel/overloads/@struct/plus.m

- Signature: `str3=plus(str1,str2)`

## Purpose

Add corresponding fields of two structures, recursively adding nested structures.

## Parameters / inputs

- `str1`, `str2` - structures with the same field names and count at each nested level.

## Outputs

- `str3` - structure containing the elementwise result of adding corresponding fields.

## Implementation structure

The function checks that both operands are structures and have matching field topology. It then applies `+` recursively to corresponding fields; a topology mismatch raises an error.
