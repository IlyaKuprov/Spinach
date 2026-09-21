# kernel/plotting/kcolourbar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/kcolourbar.m`
- Signature: `kcolourbar(x)`
- Total lines: 58

## Purpose

House style settings for Matlab figures; a product of much experience with academic publication aesthetics. Syntax: kcolourbar(x)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- x -a character string

## Outputs

- creates or updates the colour bar
- in the current axis system

## Implementation structure

- House style settings for Matlab figures; a product of much
- experience with academic publication aesthetics. Syntax:
- kcolourbar(x)
- x -a character string
- creates or updates the colour bar
- in the current axis system
- Default is empty string
- Check consistency
- Ticks to LaTeX
- Label to LaTex
- Consistency enforcement
- I've often remarked that identity politics is the product of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `colorbar()`, `ischar()`.
