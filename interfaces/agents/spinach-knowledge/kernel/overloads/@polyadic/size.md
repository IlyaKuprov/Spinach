# kernel/overloads/@polyadic/size.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/size.m`
- Signature: `varargout=size(p,dim)`
- Total lines: 81

## Purpose

Returns the size of the matrix represented by the polyadic. Syntax: answer=size(p,dim)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- p -a polyadic object
- dim -dimension whose size is required

## Outputs

- answer -a vector with one or two elements

## Implementation structure

- Returns the size of the matrix represented by the polyadic. Syntax:
- answer=size(p,dim)
- p -a polyadic object
- dim -dimension whose size is required
- answer -a vector with one or two elements
- Check consistency
- Get row dimension
- The leftmost matrix in the prefix
- The cores of the polyadic
- Get column dimension
- The rightmost matrix in the suffix
- Compose the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cellfun()`, `elseif()`, `isscalar()`, `ismember()`.
