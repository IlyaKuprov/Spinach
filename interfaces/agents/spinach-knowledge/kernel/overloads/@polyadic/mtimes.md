# kernel/overloads/@polyadic/mtimes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/overloads/@polyadic/mtimes.m`
- Signature: `C=mtimes(A,B)`
- Total lines: 171

## Purpose

Performs multiplications involving polyadics. Syntax: C=mtimes(A,B)

## Physical / mathematical content

- Polyadic tensor-product linear algebra. The emphasis is compressed operator representation, deferred algebra, and efficient Kronecker-structured manipulations.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: When A is a number; implemented by `if ~isa(A,'polyadic')&&isnumeric(A)&&isscalar(A)`.
- Lines 22-23: Multiply smallest cores in the B buffer; implemented by `for n=1:numel(B.cores)`.
- Lines 31-32: When A is a sparse matrix; implemented by `if ~isa(A,'polyadic')&&isnumeric(A)&&issparse(A)`.
- Lines 34-35: Attach as a prefix to B; implemented by `B.prefix=[{A} B.prefix]; C=simplify(B); return`.
- Lines 39-40: When A is a full matrix; implemented by `if ~isa(A,'polyadic')&&isnumeric(A)`.
- Lines 42-44: Issue a recursive call; implemented by `C=ctranspose(ctranspose(B)* ctranspose(A)); return`.
- Lines 48-49: When B is a number; implemented by `if ~isa(B,'polyadic')&&isnumeric(B)&&isscalar(B)`.
- Lines 51-52: Multiply smallest cores in the A buffer; implemented by `for n=1:numel(A.cores)`.
- Lines 60-61: When B is a sparse matrix; implemented by `if ~isa(B,'polyadic')&&isnumeric(B)&&issparse(B)`.
- Lines 63-64: Attach as a suffix to A; implemented by `A.suffix=[A.suffix {B}]; C=simplify(A); return`.
- Lines 68-69: When B is a full matrix; implemented by `if ~isa(B,'polyadic')&&isnumeric(B)`.
- Lines 71-72: Multiply by suffixes; implemented by `for n=numel(A.suffix):-1:1`.
- Lines 77-78: Preallocate the core product result; implemented by `core_rows=prod(cellfun(@(x)size(x,1),A.cores{1}))`.
- Lines 81-82: Multiply by cores; implemented by `for n=1:numel(A.cores)`.
- Lines 86-87: Multiply by prefixes; implemented by `for n=numel(A.prefix):-1:1`.
- Lines 94-95: When both are polyadic; implemented by `if isa(A,'polyadic')&&isa(B,'polyadic')`.
- Lines 97-98: Merge threshold and green light flag; implemented by `merge_thresh=1024; can_proceed=true()`.
- Lines 100-101: Does anything get in the way?; implemented by `can_proceed=isempty(A.suffix)&&can_proceed`.

### Control flow inferred from the code

- Line 20: conditional branch on `~isa(A,'polyadic')&&isnumeric(A)&&isscalar(A)`.
- Line 23: `for` loop over `n=1:numel(B.cores)`.
- Line 32: conditional branch on `~isa(A,'polyadic')&&isnumeric(A)&&issparse(A)`.
- Line 40: conditional branch on `~isa(A,'polyadic')&&isnumeric(A)`.
- Line 49: conditional branch on `~isa(B,'polyadic')&&isnumeric(B)&&isscalar(B)`.
- Line 52: `for` loop over `n=1:numel(A.cores)`.
- Line 61: conditional branch on `~isa(B,'polyadic')&&isnumeric(B)&&issparse(B)`.
- Line 69: conditional branch on `~isa(B,'polyadic')&&isnumeric(B)`.
- Line 72: `for` loop over `n=numel(A.suffix):-1:1`.
- Line 82: `for` loop over `n=1:numel(A.cores)`.
- Line 87: `for` loop over `n=numel(A.prefix):-1:1`.
- Line 95: conditional branch on `isa(A,'polyadic')&&isa(B,'polyadic')`.
- Line 111: conditional branch on `can_proceed`.
- Line 114: `for` loop over `n=1:numel(A.cores{1})`.

### Key state/data transformations

- Lines 24: computes `[~,smallest_core]` using `[~,smallest_core]=min(cellfun(@numel,B.cores{n}))`.
- Lines 25: computes `B.cores{n}{smallest_core}` using `B.cores{n}{smallest_core}=A*B.cores{n}{smallest_core}`.
- Lines 27: computes `C` using `C=simplify(B); return`.
- Lines 35: computes `B.prefix` using `B.prefix=[{A} B.prefix]; C=simplify(B); return`.
- Lines 54: computes `A.cores{n}{smallest_core}` using `A.cores{n}{smallest_core}=A.cores{n}{smallest_core}*B`.
- Lines 64: computes `A.suffix` using `A.suffix=[A.suffix {B}]; C=simplify(A); return`.
- Lines 73: computes `B` using `B=A.suffix{n}*B`.
- Lines 78: computes `core_rows` using `core_rows=prod(cellfun(@(x)size(x,1),A.cores{1}))`.
- Lines 98: computes `merge_thresh` using `merge_thresh=1024; can_proceed=true()`.
- Lines 101: computes `can_proceed` using `can_proceed=isempty(A.suffix)&&can_proceed`.
- Lines 142: computes `C.cores{1}{n}` using `C.cores{1}{n}=A.cores{1}{n}*B.cores{1}{n}`.

## Parameters / inputs

- A,B -a polyadic or a numerical array

## Outputs

- C -a polyadic or a numerical array

## Implementation structure

- Performs multiplications involving polyadics. Syntax:
- C=mtimes(A,B)
- A,B -a polyadic or a numerical array
- C -a polyadic or a numerical array
- When A is a number
- Multiply smallest cores in the B buffer
- When A is a sparse matrix
- Attach as a prefix to B
- When A is a full matrix
- Issue a recursive call
- When B is a number
- Multiply smallest cores in the A buffer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isscalar()`, `cellfun()`, `simplify()`, `issparse()`, `ctranspose()`, `kronm()`, `true()`, `polyadic()`, `suffix()`.
