# kernel/overloads/@polyadic/inflate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/inflate.m`
- Signature: `answer=inflate(p)`
- Total lines: 96

## Purpose

Converts a polyadic representation of a matrix into a sparse mat- rix. Syntax: answer=inflate(p) The function opens up all the Kronecker products while preserving the sparse type if some cores are sparse.

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Process nested polyadics; implemented by `for n=1:numel(p.cores)`.
- Lines 43-44: Find the core dimensions; implemented by `nrows=prod(cellfun(@(x)size(x,1),p.cores{1}))`.
- Lines 47-48: Get index arrays going; implemented by `rows=cell(numel(p.cores),1)`.
- Lines 52-53: Loop over the sum; implemented by `for n=1:numel(p.cores)`.
- Lines 55-56: Compute the polyadic; implemented by `term=p.cores{n}{1}`.
- Lines 64-65: Merge index arrays; implemented by `rows=cell2mat(rows)`.
- Lines 69-70: Build a sparse matrix; implemented by `if isempty(vals)`.
- Lines 76-77: Prevent GPU memory leaks; implemented by `clear('rows','cols','vals','term')`.
- Lines 79-80: Multiply by prefixes; implemented by `for n=numel(p.prefix):-1:1`.
- Lines 84-85: Multiply by suffixes; implemented by `for n=1:numel(p.suffix)`.

### Control flow inferred from the code

- Line 25: `for` loop over `n=1:numel(p.cores)`.
- Line 26: `for` loop over `k=1:numel(p.cores{n})`.
- Line 27: conditional branch on `isa(p.cores{n}{k},'polyadic')`.
- Line 32: `for` loop over `n=1:numel(p.prefix)`.
- Line 33: conditional branch on `isa(p.prefix{n},'polyadic')`.
- Line 37: `for` loop over `n=1:numel(p.suffix)`.
- Line 38: conditional branch on `isa(p.suffix{n},'polyadic')`.
- Line 53: `for` loop over `n=1:numel(p.cores)`.
- Line 57: `for` loop over `k=2:numel(p.cores{n})`.
- Line 70: conditional branch on `isempty(vals)`.
- Line 80: `for` loop over `n=numel(p.prefix):-1:1`.
- Line 85: `for` loop over `n=1:numel(p.suffix)`.

### Key state/data transformations

- Lines 28: computes `p.cores{n}{k}` using `p.cores{n}{k}=inflate(p.cores{n}{k})`.
- Lines 34: computes `p.prefix{n}` using `p.prefix{n}=inflate(p.prefix{n})`.
- Lines 39: computes `p.suffix{n}` using `p.suffix{n}=inflate(p.suffix{n})`.
- Lines 44: computes `nrows` using `nrows=prod(cellfun(@(x)size(x,1),p.cores{1}))`.
- Lines 45: computes `ncols` using `ncols=prod(cellfun(@(x)size(x,2),p.cores{1}))`.
- Lines 48: computes `rows` using `rows=cell(numel(p.cores),1)`.
- Lines 49: computes `cols` using `cols=cell(numel(p.cores),1)`.
- Lines 50: computes `vals` using `vals=cell(numel(p.cores),1)`.
- Lines 56: computes `term` using `term=p.cores{n}{1}`.
- Lines 60: computes `[rows{n},cols{n},vals{n}]` using `[rows{n},cols{n},vals{n}]=find(term)`.
- Lines 71: computes `answer` using `answer=sparse([],[],[],nrows,ncols)`.

## Parameters / inputs

- p -a polyadic object

## Outputs

- answer -a sparse matrix, except if prefixes or suffixes
- are full (in that case, a full matrix)

## Implementation structure

- Converts a polyadic representation of a matrix into a sparse mat-
- rix. Syntax:
- answer=inflate(p)
- The function opens up all the Kronecker products while preserving
- the sparse type if some cores are sparse.
- p -a polyadic object
- answer -a sparse matrix, except if prefixes or suffixes
- are full (in that case, a full matrix)
- Process nested polyadics
- Find the core dimensions
- Get index arrays going
- Loop over the sum

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`, `cell2mat()`, `clear()`.
