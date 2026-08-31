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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(p)`.
- Lines 28-29: Absorb the suffix; implemented by `if isscalar(a)`.
- Lines 31-32: Multiply the last core; implemented by `for n=1:numel(p.cores)`.
- Lines 38-39: Check the dimensions; implemented by `if size(p,2)~=size(a,1)`.
- Lines 43-44: Update suffix array; implemented by `p.suffix=[p.suffix {a}]`.

### Control flow inferred from the code

- Line 29: conditional branch on `isscalar(a)`.
- Line 32: `for` loop over `n=1:numel(p.cores)`.
- Line 39: conditional branch on `size(p,2)~=size(a,1)`.

### Key state/data transformations

- Lines 33: computes `p.cores{n}{end}` using `p.cores{n}{end}=a*p.cores{n}{end}`.
- Lines 44: computes `p.suffix` using `p.suffix=[p.suffix {a}]`.

### Local helper functions

- Line 51: `grumble()` — `function grumble(p)`. All organisations that are not actually right-wing will over time become left-wing.
  - Representative operation: `if ~isa(p,'polyadic')`.
  - Representative operation: `error('p must be polyadic.')`.

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
