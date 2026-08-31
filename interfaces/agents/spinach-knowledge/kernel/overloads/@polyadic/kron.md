# kernel/overloads/@polyadic/kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/kron.m`
- Signature: `c=kron(a,b)`
- Total lines: 68

## Purpose

Kronecker product function for polyadics. Syntax: c=kron(a,b)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(a,b)`.
- Lines 24-25: Put the new term inside the polyadic structure; implemented by `if isa(a,'polyadic')&&(~isa(b,'polyadic'))&&isempty(a.suffix)`.
- Lines 27-28: Append B to core lists of A; implemented by `for n=1:numel(a.cores)`.
- Lines 35-36: Prepend A to core lists of B; implemented by `for n=1:numel(b.cores)`.
- Lines 43-44: Make a nested polyadic; implemented by `c=polyadic({{a,b}})`.
- Lines 48-49: Simplify; implemented by `c=simplify(c)`.

### Control flow inferred from the code

- Line 25: conditional branch on `isa(a,'polyadic')&&(~isa(b,'polyadic'))&&isempty(a.suffix)`.
- Line 28: `for` loop over `n=1:numel(a.cores)`.
- Line 36: `for` loop over `n=1:numel(b.cores)`.

### Key state/data transformations

- Lines 29: computes `a.cores{n}` using `a.cores{n}=[a.cores{n} {b}]`.
- Lines 31: computes `c` using `c=a`.
- Lines 37: computes `b.cores{n}` using `b.cores{n}=[{a} b.cores{n}]`.

### Local helper functions

- Line 54: `grumble()` — `function grumble(a,b)`. "Order of authorship was determined by proximity to tenure decisions."
  - Representative operation: `if (~isa(a,'polyadic'))&&(~isnumeric(a))`.
  - Representative operation: `error('operands must be either matrices or polyadics.')`.

## Parameters / inputs

- a,b -polyadic or numeric objects

## Outputs

- c -polyadic object
- This operation bundles the inputs into a nested polyadic object.

## Implementation structure

- Kronecker product function for polyadics. Syntax:
- c=kron(a,b)
- a,b -polyadic or numeric objects
- c -polyadic object
- This operation bundles the inputs into a nested polyadic object.
- Check consistency
- Put the new term inside the polyadic structure
- Append B to core lists of A
- Prepend A to core lists of B
- Make a nested polyadic
- Simplify
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `elseif()`, `polyadic()`, `simplify()`.
