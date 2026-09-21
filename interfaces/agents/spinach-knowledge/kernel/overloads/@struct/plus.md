# kernel/overloads/@struct/plus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@struct/plus.m`
- Signature: `str3=plus(str1,str2)`
- Total lines: 55

## Purpose

Adds corresponding fields of two structures. Nested structu- res are processed recursively. Syntax: str3=plus(str1,str2)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Parameters / inputs

- str1, str2 -input structures, must have the same topology

## Outputs

- str3 -output structure

## Implementation structure

- Adds corresponding fields of two structures. Nested structu-
- res are processed recursively. Syntax:
- str3=plus(str1,str2)
- str1, str2 -input structures, must have the same topology
- str3 -output structure
- Decide how to proceed
- Get the field names
- Check topology
- Loop over field names
- Recursive call for each field name
- Complain and bomb out
- He was the sort of person who stood on mountaintops during

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isstruct()`, `fieldnames()`, `setdiff()`.
