# kernel/overloads/@polyadic/validate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/validate.m`
- Signature: `validate(p)`
- Total lines: 105

## Purpose

Checks the internal structure of a polyadic object and throws an error if the object does not meet expectations. Syntax: validate(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Check the type; implemented by `if ~isa(p,'polyadic'), error('this object is not a polyadic.'); end`.
- Lines 19-20: Check core array, level 1; implemented by `if ~iscell(p.cores), error('p.cores must be a cell array.'); end`.
- Lines 22-23: Check prefix and suffix arrays, level 1; implemented by `if ~iscell(p.prefix), error('p.prefix must be a cell array.'); end`.
- Lines 26-27: Check core array, level 2; implemented by `for n=1:numel(p.cores)`.
- Lines 39-40: Check prefix and suffix arrays, level 2; implemented by `for n=1:numel(p.prefix)`.
- Lines 53-54: Check core dimensions; implemented by `core_dims=zeros(numel(p.cores),2)`.
- Lines 66-67: Check prefix and suffix dimensions; implemented by `if (~isempty(p.prefix))&&(size(p.prefix{end},2)~=core_dims(1,1))`.
- Lines 88-89: Check the number of terms; implemented by `if numel(p.cores)>100`.

### Control flow inferred from the code

- Line 17: conditional branch on `~isa(p,'polyadic'), error('this object is not a polyadic.'); end`.
- Line 20: conditional branch on `~iscell(p.cores), error('p.cores must be a cell array.'); end`.
- Line 23: conditional branch on `~iscell(p.prefix), error('p.prefix must be a cell array.'); end`.
- Line 24: conditional branch on `~iscell(p.suffix), error('p.suffix must be a cell array.'); end`.
- Line 27: `for` loop over `n=1:numel(p.cores)`.
- Line 28: conditional branch on `~iscell(p.cores{n})`.
- Line 31: `for` loop over `k=1:numel(p.cores{n})`.
- Line 32: conditional branch on `~isnumeric(p.cores{n}{k})`.
- Line 35: conditional branch on `isa(p.cores{n}{k},'polyadic'), validate(p.cores{n}{k}); end`.
- Line 40: `for` loop over `n=1:numel(p.prefix)`.
- Line 41: conditional branch on `~isnumeric(p.prefix{n})`.
- Line 44: conditional branch on `isa(p.prefix{n},'polyadic'), validate(p.prefix{n}); end`.
- Line 46: `for` loop over `n=1:numel(p.suffix)`.
- Line 47: conditional branch on `~isnumeric(p.suffix{n})`.

### Key state/data transformations

- Lines 54: computes `core_dims` using `core_dims=zeros(numel(p.cores),2)`.
- Lines 56: computes `nrows` using `nrows=cellfun(@(x)size(x,1),p.cores{n})`.
- Lines 57: computes `ncols` using `ncols=cellfun(@(x)size(x,2),p.cores{n})`.
- Lines 58: computes `core_dims(n,1)` using `core_dims(n,1)=prod(nrows(:))`.
- Lines 59: computes `core_dims(n,2)` using `core_dims(n,2)=prod(ncols(:))`.

## Parameters / inputs

- p -a polyadic object

## Implementation structure

- Checks the internal structure of a polyadic object and throws an
- error if the object does not meet expectations. Syntax:
- validate(p)
- p -a polyadic object
- Check the type
- Check core array, level 1
- Check prefix and suffix arrays, level 1
- Check core array, level 2
- Check prefix and suffix arrays, level 2
- Check core dimensions
- Check prefix and suffix dimensions
- Check the number of terms

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `iscell()`, `cellfun()`, `core_dims()`, `nrows()`, `ncols()`, `all()`.
