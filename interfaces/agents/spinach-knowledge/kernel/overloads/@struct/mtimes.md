# kernel/overloads/@struct/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@struct/mtimes.m`
- Signature: `str_out=mtimes(M,str_in)`
- Total lines: 56

## Purpose

Multiplies all entries of a structure by a user-specified mat- rix. Nested structures are processed recursively. Syntax: str_out=mtimes(M,str_in)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fieldnames()`, `isstruct()`.
