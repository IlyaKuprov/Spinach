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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Default is empty string; implemented by `if ~exist('x','var'), x=''; end`.
- Lines 24-25: Check consistency; implemented by `grumble(x)`.
- Lines 27-29: Ticks to LaTeX; implemented by `cb=colorbar('TickLabelInterpreter','latex', 'FontSize',12)`.
- Lines 31-32: Label to LaTex; implemented by `cb.Label.Interpreter='latex'`.

### Control flow inferred from the code

- Line 22: conditional branch on `~exist('x','var'), x=''; end`.

### Key state/data transformations

- Lines 28-29: computes `cb` using `cb=colorbar('TickLabelInterpreter','latex', 'FontSize',12)`.
- Lines 32: computes `cb.Label.Interpreter` using `cb.Label.Interpreter='latex'`.
- Lines 33: computes `cb.Label.FontSize` using `cb.Label.FontSize=13`.
- Lines 34: computes `cb.Label.String` using `cb.Label.String=x`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(x)`. I've often remarked that identity politics is the product of prosperity. The movement's nitpicking about undetectable-to-
  - Representative operation: `if ~ischar(x)`.
  - Representative operation: `error('x must be a character string')`.

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
