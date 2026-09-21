# kernel/overloads/@polyadic/suffix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/suffix.m`
- Signature: `p=suffix(p,a)`
- Total lines: 61

## Purpose

Adds suffix matrices to a polyadic. Anything the polyadic multiplies will first be multiplied by the suffix matri- ces. Syntax: p=suffix(p,a)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- p -polyadic object
- a -suffix matrix

## Outputs

- p -polyadic object
- Note: a suffix can be a polyadic itself.

## Implementation structure

- Adds suffix matrices to a polyadic. Anything the polyadic
- multiplies will first be multiplied by the suffix matri-
- ces. Syntax:
- p=suffix(p,a)
- p - polyadic object
- a - suffix matrix
- Note: a suffix can be a polyadic itself.
- Check consistency
- Absorb the suffix
- Multiply the last core
- Check the dimensions
- Update suffix array

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
