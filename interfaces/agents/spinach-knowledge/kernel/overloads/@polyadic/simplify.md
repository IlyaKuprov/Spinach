# kernel/overloads/@polyadic/simplify.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/simplify.m`
- Signature: `p=simplify(p)`
- Total lines: 190

## Purpose

Simplifies the structure of the polyadic object by reordering buffers, dropping inconsequential terms, and flattening nested polyadics where possible. Syntax: p=simplify(p)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Check consistency; implemented by `grumble(p)`.
- Lines 24-25: Get size information; implemented by `[nrows,ncols]=size(p)`.
- Lines 27-28: Flush trivial polyadics into all-zero sparse matrices; implemented by `if isempty(p.cores), p=spalloc(nrows,ncols,0); return; end`.
- Lines 30-31: Loop until static; implemented by `changes_made=true()`.
- Lines 34-35: Default disposition; implemented by `changes_made=false()`.
- Lines 37-38: Simplify prefixes; implemented by `for n=1:numel(p.prefix)`.
- Lines 40-41: Recursive call for polyadics; implemented by `if isa(p.prefix{n},'polyadic')`.
- Lines 45-46: A zero prefix flushes into all-zero sparse matrix; implemented by `if nnz(p.prefix{n})==0, p=spalloc(nrows,ncols,0); return; end`.
- Lines 48-49: Drop unit prefixes; implemented by `if iseye(p.prefix{n})`.
- Lines 54-55: Absorb opium prefixes; implemented by `if isa(p.prefix{n},'opium')`.
- Lines 60-61: Absorb scalar prefixes; implemented by `if isnumeric(p.prefix{n})&&isscalar(p.prefix{n})`.
- Lines 69-70: Simplify suffixes; implemented by `for n=1:numel(p.suffix)`.
- Lines 72-73: Recursive call for polyadics; implemented by `if isa(p.suffix{n},'polyadic')`.
- Lines 77-78: A zero suffix flushes into all-zero sparse matrix; implemented by `if nnz(p.suffix{n})==0, p=spalloc(nrows,ncols,0); return; end`.
- Lines 80-81: Drop unit suffixes; implemented by `if iseye(p.suffix{n})`.
- Lines 86-87: Absorb opium suffixes; implemented by `if isa(p.suffix{n},'opium')`.
- Lines 92-93: Absorb scalar suffixes; implemented by `if isnumeric(p.suffix{n})&&isscalar(p.suffix{n})`.
- Lines 101-102: Over sum terms; implemented by `for n=1:numel(p.cores)`.

### Control flow inferred from the code

- Line 28: conditional branch on `isempty(p.cores), p=spalloc(nrows,ncols,0); return; end`.
- Line 32: `while` loop over `changes_made`.
- Line 38: `for` loop over `n=1:numel(p.prefix)`.
- Line 41: conditional branch on `isa(p.prefix{n},'polyadic')`.
- Line 46: conditional branch on `nnz(p.prefix{n})==0, p=spalloc(nrows,ncols,0); return; end`.
- Line 49: conditional branch on `iseye(p.prefix{n})`.
- Line 55: conditional branch on `isa(p.prefix{n},'opium')`.
- Line 61: conditional branch on `isnumeric(p.prefix{n})&&isscalar(p.prefix{n})`.
- Line 70: `for` loop over `n=1:numel(p.suffix)`.
- Line 73: conditional branch on `isa(p.suffix{n},'polyadic')`.
- Line 78: conditional branch on `nnz(p.suffix{n})==0, p=spalloc(nrows,ncols,0); return; end`.
- Line 81: conditional branch on `iseye(p.suffix{n})`.
- Line 87: conditional branch on `isa(p.suffix{n},'opium')`.
- Line 93: conditional branch on `isnumeric(p.suffix{n})&&isscalar(p.suffix{n})`.

### Key state/data transformations

- Lines 25: computes `[nrows,ncols]` using `[nrows,ncols]=size(p)`.
- Lines 31: computes `changes_made` using `changes_made=true()`.
- Lines 42: computes `p.prefix{n}` using `p.prefix{n}=simplify(p.prefix{n})`.
- Lines 56: computes `coeff` using `coeff=p.prefix{n}.coeff`.
- Lines 57: computes `p.prefix(n)` using `p.prefix(n)=[]; p=coeff*p; return`.
- Lines 67: computes `p.prefix(cellfun(@isempty,p.prefix))` using `p.prefix(cellfun(@isempty,p.prefix))=[]`.
- Lines 74: computes `p.suffix{n}` using `p.suffix{n}=simplify(p.suffix{n})`.
- Lines 89: computes `p.suffix(n)` using `p.suffix(n)=[]; p=coeff*p; return`.
- Lines 99: computes `p.suffix(cellfun(@isempty,p.suffix))` using `p.suffix(cellfun(@isempty,p.suffix))=[]`.
- Lines 111: computes `p.cores{n}{k}` using `p.cores{n}{k}=simplify(p.cores{n}{k})`.
- Lines 120-122: computes `p.cores{n}` using `p.cores{n}=[p.cores{n}(1:(k-1)) p.cores{n}{k}.cores{1} p.cores{n}((k+1):end)]`.
- Lines 143: computes `p.cores(cellfun(@isempty,p.cores))` using `p.cores(cellfun(@isempty,p.cores))=[]`.
- Lines 159: computes `p.cores{n}(k+1)` using `p.cores{n}(k+1)=[]; changes_made=true(); break`.
- Lines 172: computes `p` using `p=p.cores{1}{1}`.

### Local helper functions

- Line 178: `grumble()` — `function grumble(p)`. New scientific findings do not become commonly accepted by convincing other scientists, but rather by waiting until
  - Representative operation: `if ~isa(p,'polyadic')`.
  - Representative operation: `error('p must be polyadic.')`.

## Parameters / inputs

- p -a polyadic object

## Outputs

- p -a polyadic or a numeric object

## Implementation structure

- Simplifies the structure of the polyadic object by reordering buffers,
- dropping inconsequential terms, and flattening nested polyadics where
- possible. Syntax:
- p=simplify(p)
- p -a polyadic object
- p -a polyadic or a numeric object
- Check consistency
- Get size information
- Flush trivial polyadics into all-zero sparse matrices
- Loop until static
- Default disposition
- Simplify prefixes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `true()`, `false()`, `nnz()`, `iseye()`, `isscalar()`, `cellfun()`, `opium()`.
