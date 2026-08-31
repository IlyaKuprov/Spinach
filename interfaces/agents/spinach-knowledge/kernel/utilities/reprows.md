# kernel/utilities/reprows.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/reprows.m`
- Signature: `B=reprows(A,row_nums,rep_counts)`
- Total lines: 67

## Purpose

Replicates specified rows of a matrix or cell array a specified number of times. Syntax: B=reprows(A,row_nums,rep_counts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(A,row_nums,rep_counts)`.
- Lines 28-29: Replication counts for every row; implemented by `n=size(A,1); rep_map=ones(1,n)`.
- Lines 32-33: Build row index vector; implemented by `row_idx=repelem(1:n,rep_map)`.
- Lines 35-36: Extract and replicate; implemented by `B=A(row_idx,:)`.

### Key state/data transformations

- Lines 29: computes `n` using `n=size(A,1); rep_map=ones(1,n)`.
- Lines 30: computes `rep_map(row_nums)` using `rep_map(row_nums)=rep_counts(:)`.
- Lines 33: computes `row_idx` using `row_idx=repelem(1:n,rep_map)`.
- Lines 36: computes `B` using `B=A(row_idx,:)`.

### Local helper functions

- Line 41: `grumble()` — `function grumble(A,row_nums,rep_counts)`.
  - Representative operation: `if (~isnumeric(A))&&(~iscell(A))`.
  - Representative operation: `error('A must be numeric or a cell array.')`.

## Parameters / inputs

- A -a numeric matrix or a cell array
- row_nums -vector of row indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each row to make

## Outputs

- B -same type as A

## Implementation structure

- Replicates specified rows of a matrix or cell array a
- specified number of times. Syntax:
- B=reprows(A,row_nums,rep_counts)
- A -a numeric matrix or a cell array
- row_nums -vector of row indices to replicate
- rep_counts -vector of positive integers specifying
- how many copies of each row to make
- B -same type as A
- Check consistency
- Replication counts for every row
- Build row index vector
- Extract and replicate

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rep_map()`, `rep_counts()`, `repelem()`, `iscell()`, `isvector()`, `any()`.
