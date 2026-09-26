# kernel/overloads/@struct/mtimes.m

- Signature: `str_out=mtimes(M,str_in)`

## Purpose

Multiplies all entries of a structure by a user-specified mat- rix. Nested structures are processed recursively. Syntax: str_out=mtimes(M,str_in)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- M -any numeric object (scalar, matrix, etc.)
- str_in -a structure with numeric subfields

## Outputs

- str_out -the resulting structure

## Implementation structure

- Multiplies all entries of a structure by a user-specified mat-
- rix. Nested structures are processed recursively. Syntax:
- str_out=mtimes(M,str_in)
- M -any numeric object (scalar, matrix, etc.)
- str_in -a structure with numeric subfields
- str_out -the resulting structure
- Check consistency
- Get the field names
- Loop over field names
- Recursive call for each field name
- Consistency enforcement
- Arthur Dent: What happens if I press this button?
