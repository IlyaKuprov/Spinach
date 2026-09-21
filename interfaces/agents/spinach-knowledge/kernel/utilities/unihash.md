# kernel/utilities/unihash.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/unihash.m`
- Signature: `A=unihash(A)`
- Total lines: 55

## Purpose

Hash table based stable duplicate row eliminator, for use with large sparse matrices where Matlab's unique(...,'rows') is too slow. Syntax: A=unihash(A)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -a large and sparse matrix

## Outputs

- A -same matrix with duplicate
- rows deleted, keeping the
- first occurrence of each

## Implementation structure

- Hash table based stable duplicate row eliminator,
- for use with large sparse matrices where Matlab's
- unique(...,'rows') is too slow. Syntax:
- A=unihash(A)
- A -a large and sparse matrix
- A -same matrix with duplicate
- rows deleted, keeping the
- first occurrence of each
- Check consistency
- Build an MD5 hash table
- Redundant row index using a hash table
- Elimination

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hash_table()`, `md5_hash()`, `ismatrix()`.
