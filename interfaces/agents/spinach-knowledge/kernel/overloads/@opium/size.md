# kernel/overloads/@opium/size.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@opium/size.m`
- Signature: `varargout=size(op,dim)`
- Total lines: 52

## Purpose

The size of the matrix represented by the OPIUM. Syntax: answer=size(op,dim)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- op -an opium object
- dim -optional, dimension whose
- size is required

## Outputs

- answer -a vector with one or two elements

## Implementation structure

- The size of the matrix represented by the OPIUM. Syntax:
- answer=size(op,dim)
- op -an opium object
- dim -optional, dimension whose
- size is required
- answer -a vector with one or two elements
- Check consistency
- Compose the answer
- Consistency enforcement
- Treason doth never prosper: what's the reason?
- Why, if it prosper, none dare call it treason.
- John Harrington

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `elseif()`, `isscalar()`, `ismember()`.
