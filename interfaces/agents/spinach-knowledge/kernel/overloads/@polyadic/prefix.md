# kernel/overloads/@polyadic/prefix.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/prefix.m`
- Signature: `p=prefix(a,p)`
- Total lines: 62

## Purpose

Adds prefix matrices to a polyadic. Anything the polyadic multiplies will subsequently be multiplied by the prefix matrices. Syntax: p=prefix(a,p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- a -prefix matrix
- p -polyadic object

## Outputs

- p -polyadic object
- Note: a prefix can be a polyadic itself.

## Implementation structure

- Adds prefix matrices to a polyadic. Anything the polyadic
- multiplies will subsequently be multiplied by the prefix
- matrices. Syntax:
- p=prefix(a,p)
- a - prefix matrix
- p - polyadic object
- Note: a prefix can be a polyadic itself.
- Check consistency
- Absorb the prefix
- Multiply the first core
- Check the dimensions
- Update prefix array

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
