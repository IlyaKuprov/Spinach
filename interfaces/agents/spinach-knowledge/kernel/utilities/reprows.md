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
