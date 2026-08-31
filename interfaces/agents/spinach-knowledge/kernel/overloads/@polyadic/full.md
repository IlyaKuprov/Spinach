# kernel/overloads/@polyadic/full.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/full.m`
- Signature: `answer=full(p)`
- Total lines: 79

## Purpose

Converts a polyadic representation of a matrix into a full mat- rix. Syntax: answer=full(p) The function opens up all the Kronecker products and uses full arithmetic throughout even if some cores are sparse.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Process nested polyadics; implemented by `for n=1:numel(p.cores)`.
- Lines 42-43: Find the core dimensions; implemented by `nrows=prod(cellfun(@(x)size(x,1),p.cores{1}))`.
- Lines 46-47: Preallocate the answer; implemented by `answer=zeros(nrows,ncols)`.
- Lines 49-50: Loop over the sum; implemented by `for n=1:numel(p.cores)`.
- Lines 52-53: Compute the polyadic; implemented by `term=full(p.cores{n}{1})`.
- Lines 61-62: Multiply by prefixes; implemented by `for n=numel(p.prefix):-1:1`.
- Lines 66-67: Multiply by suffixes; implemented by `for n=1:numel(p.suffix)`.
- Lines 71-72: Make sure the result is full; implemented by `answer=full(answer)`.

### Control flow inferred from the code

- Line 24: `for` loop over `n=1:numel(p.cores)`.
- Line 25: `for` loop over `k=1:numel(p.cores{n})`.
- Line 26: conditional branch on `isa(p.cores{n}{k},'polyadic')`.
- Line 31: `for` loop over `n=1:numel(p.prefix)`.
- Line 32: conditional branch on `isa(p.prefix{n},'polyadic')`.
- Line 36: `for` loop over `n=1:numel(p.suffix)`.
- Line 37: conditional branch on `isa(p.suffix{n},'polyadic')`.
- Line 50: `for` loop over `n=1:numel(p.cores)`.
- Line 54: `for` loop over `k=2:numel(p.cores{n})`.
- Line 62: `for` loop over `n=numel(p.prefix):-1:1`.
- Line 67: `for` loop over `n=1:numel(p.suffix)`.

### Key state/data transformations

- Lines 27: computes `p.cores{n}{k}` using `p.cores{n}{k}=full(p.cores{n}{k})`.
- Lines 33: computes `p.prefix{n}` using `p.prefix{n}=full(p.prefix{n})`.
- Lines 38: computes `p.suffix{n}` using `p.suffix{n}=full(p.suffix{n})`.
- Lines 43: computes `nrows` using `nrows=prod(cellfun(@(x)size(x,1),p.cores{1}))`.
- Lines 44: computes `ncols` using `ncols=prod(cellfun(@(x)size(x,2),p.cores{1}))`.
- Lines 47: computes `answer` using `answer=zeros(nrows,ncols)`.
- Lines 53: computes `term` using `term=full(p.cores{n}{1})`.

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -a full matrix

## Implementation structure

- Converts a polyadic representation of a matrix into a full mat-
- rix. Syntax:
- answer=full(p)
- The function opens up all the Kronecker products and uses full
- arithmetic throughout even if some cores are sparse.
- p -a polyadic object
- answer -a full matrix
- Process nested polyadics
- Find the core dimensions
- Preallocate the answer
- Loop over the sum
- Compute the polyadic

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`.
