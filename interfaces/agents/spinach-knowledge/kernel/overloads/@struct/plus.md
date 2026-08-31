# kernel/overloads/@struct/plus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@struct/plus.m`
- Signature: `str3=plus(str1,str2)`
- Total lines: 55

## Purpose

Adds corresponding fields of two structures. Nested structu- res are processed recursively. Syntax: str3=plus(str1,str2)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Decide how to proceed; implemented by `if isstruct(str1)&&isstruct(str2)`.
- Lines 24-25: Get the field names; implemented by `fnames1=fieldnames(str1); fnames2=fieldnames(str2)`.
- Lines 27-29: Check topology; implemented by `if (numel(fnames1)~=numel(fnames2))|| (~isempty(setdiff(fnames1,fnames2)))`.
- Lines 33-34: Loop over field names; implemented by `for n=1:length(fnames1)`.
- Lines 36-37: Recursive call for each field name; implemented by `str3.(fnames1{n})=str1.(fnames1{n})+str2.(fnames1{n})`.
- Lines 43-44: Complain and bomb out; implemented by `error('structure topology mismatch.')`.

### Control flow inferred from the code

- Line 22: conditional branch on `isstruct(str1)&&isstruct(str2)`.
- Line 28: conditional branch on `(numel(fnames1)~=numel(fnames2))||`.
- Line 34: `for` loop over `n=1:length(fnames1)`.

### Key state/data transformations

- Lines 25: computes `fnames1` using `fnames1=fieldnames(str1); fnames2=fieldnames(str2)`.
- Lines 37: computes `str3.(fnames1{n})` using `str3.(fnames1{n})=str1.(fnames1{n})+str2.(fnames1{n})`.

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
