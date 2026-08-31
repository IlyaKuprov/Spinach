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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `if nargin==2, grumble(dim); end`.
- Lines 24-25: Get row dimension; implemented by `if ~isempty(p.prefix)`.
- Lines 27-28: The leftmost matrix in the prefix; implemented by `nrows=size(p.prefix{1},1)`.
- Lines 32-33: The cores of the polyadic; implemented by `nrows=prod(cellfun(@(x)size(x,1),p.cores{1}))`.
- Lines 37-38: Get column dimension; implemented by `if ~isempty(p.suffix)`.
- Lines 40-41: The rightmost matrix in the suffix; implemented by `ncols=size(p.suffix{end},2)`.
- Lines 45-46: The cores of the polyadic; implemented by `ncols=prod(cellfun(@(x)size(x,2),p.cores{1}))`.
- Lines 50-51: Compose the answer; implemented by `if (nargin==1)&&(nargout<=1)`.

### Control flow inferred from the code

- Line 22: conditional branch on `nargin==2, grumble(dim); end`.
- Line 25: conditional branch on `~isempty(p.prefix)`.
- Line 38: conditional branch on `~isempty(p.suffix)`.
- Line 51: conditional branch on `(nargin==1)&&(nargout<=1)`.

### Key state/data transformations

- Lines 28: computes `nrows` using `nrows=size(p.prefix{1},1)`.
- Lines 41: computes `ncols` using `ncols=size(p.suffix{end},2)`.
- Lines 52: computes `varargout{1}` using `varargout{1}=[nrows ncols]`.
- Lines 55: computes `varargout{2}` using `varargout{2}=ncols`.

### Local helper functions

- Line 67: `grumble()` — `function grumble(dim)`. What is true, just and beautiful is not determined by popular vote. The masses everywhere are ignorant, short-sighted, moti-
  - Representative operation: `if (~isscalar(dim))||(~ismember(dim,[1 2]))`.
  - Representative operation: `error('for a polyadic object, dim must be 1 or 2')`.

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
