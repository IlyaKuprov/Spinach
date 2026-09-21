# kernel/utilities/polinfo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/polinfo.m`
- Signature: `polinfo(p,level,label)`
- Total lines: 111

## Purpose

Draws an ASCII diagram of a polyadic object. Syntax: polinfo(p)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- p -polyadic object

## Outputs

- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.

## Implementation structure

- Draws an ASCII diagram of a polyadic object. Syntax:
- polinfo(p)
- p - polyadic object
- an ASCII diagram to the console
- Note: polyadic objects can be huge, the code below
- avoids making memory copies.
- Default settings
- Check consistency
- Set indentation level
- Get the dimensions
- Print label and size
- Print prefixes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `isscalar()`, `ischar()`.
