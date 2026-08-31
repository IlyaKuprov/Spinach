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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(p)`.
- Lines 28-29: Absorb the prefix; implemented by `if isscalar(a)`.
- Lines 31-32: Multiply the first core; implemented by `for n=1:numel(p.cores)`.
- Lines 38-39: Check the dimensions; implemented by `if size(a,2)~=size(p,1)`.
- Lines 43-44: Update prefix array; implemented by `p.prefix=[{a} p.prefix]`.

### Control flow inferred from the code

- Line 29: conditional branch on `isscalar(a)`.
- Line 32: `for` loop over `n=1:numel(p.cores)`.
- Line 39: conditional branch on `size(a,2)~=size(p,1)`.

### Key state/data transformations

- Lines 33: computes `p.cores{n}{1}` using `p.cores{n}{1}=a*p.cores{n}{1}`.
- Lines 44: computes `p.prefix` using `p.prefix=[{a} p.prefix]`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(p)`. It is difficult to get a man to understand something when his salary depends on his not
  - Representative operation: `if ~isa(p,'polyadic')`.
  - Representative operation: `error('p must be polyadic.')`.

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
