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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Check consistency; implemented by `grumble(M,str_in)`.
- Lines 26-27: Get the field names; implemented by `fnames=fieldnames(str_in)`.
- Lines 29-30: Loop over field names; implemented by `for n=1:numel(fnames)`.
- Lines 32-33: Recursive call for each field name; implemented by `str_out.(fnames{n})=M*str_in.(fnames{n})`.

### Control flow inferred from the code

- Line 30: `for` loop over `n=1:numel(fnames)`.

### Key state/data transformations

- Lines 27: computes `fnames` using `fnames=fieldnames(str_in)`.
- Lines 33: computes `str_out.(fnames{n})` using `str_out.(fnames{n})=M*str_in.(fnames{n})`.

### Local helper functions

- Line 40: `grumble()` — `function grumble(M,str_in)`. Arthur Dent: What happens if I press this button? Ford Prefect: I wouldn't --
  - Representative operation: `if ~isnumeric(M)`.
  - Representative operation: `error('the first argument must be numeric.')`.

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
