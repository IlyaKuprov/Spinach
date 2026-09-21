# kernel/plotting/scale_figure.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/scale_figure.m`
- Signature: `scale_figure(by)`
- Total lines: 55

## Purpose

Scales the current figure from the default size by the factors provided by the user. Syntax: scale_figure(by)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- by -two-element row vector of scaling factors,
- format: [width height]

## Implementation structure

- Scales the current figure from the default size by the factors
- provided by the user. Syntax:
- scale_figure(by)
- by - two-element row vector of scaling factors,
- format: [width height]
- Check consistency
- Get figure location
- Get figure centroid
- Get default figure size
- Scale the figure
- Update figure parameters
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `get()`, `loc()`, `set()`, `isrow()`, `any()`.
