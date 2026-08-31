# kernel/utilities/repcols.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/repcols.m`
- Signature: `B=repcols(A,col_nums,rep_counts)`
- Total lines: 66

## Purpose

Replicates specified columns of a matrix or cell array a specified number of times. Syntax: B=repcols(A,col_nums,rep_counts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(A,col_nums,rep_counts)`.
- Lines 28-29: Replication counts for every column; implemented by `n=size(A,2); rep_map=ones(1,n)`.
- Lines 32-33: Build column index vector; implemented by `col_idx=repelem(1:n,rep_map)`.
- Lines 35-36: Extract and replicate; implemented by `B=A(:,col_idx)`.

### Key state/data transformations

- Lines 29: computes `n` using `n=size(A,2); rep_map=ones(1,n)`.
- Lines 30: computes `rep_map(col_nums)` using `rep_map(col_nums)=rep_counts(:)`.
- Lines 33: computes `col_idx` using `col_idx=repelem(1:n,rep_map)`.
- Lines 36: computes `B` using `B=A(:,col_idx)`.

### Local helper functions

- Line 41: `grumble()` — `function grumble(A,col_nums,rep_counts)`.
  - Representative operation: `if (~isnumeric(A))&&(~iscell(A))`.
  - Representative operation: `error('A must be numeric or a cell array.')`.

## Parameters / inputs

- A -a numeric matrix or a cell array
- col_nums -vector of column indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each column to make
- Output:
- B -same type as A

## Implementation structure

- Replicates specified columns of a matrix or cell array a
- specified number of times. Syntax:
- B=repcols(A,col_nums,rep_counts)
- A -a numeric matrix or a cell array
- col_nums -vector of column indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each column to make
- Output:
- B -same type as A
- Check consistency
- Replication counts for every column
- Build column index vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rep_map()`, `rep_counts()`, `repelem()`, `iscell()`, `isvector()`, `any()`.
