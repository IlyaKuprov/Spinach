# kernel/utilities/clean_up.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/clean_up.m`
- Signature: `A=clean_up(spin_system,A,nonzero_tol)`
- Total lines: 98

## Purpose

Array clean-up utility. Drops non-zero elements with magnitude below the user-specified tolerance and converts between sparse and full storage de- pending on the density of non-zeroes in the array. Syntax: A=clean_up(spin_system,A,nonzero_tol)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -a numerical array or a cell array thereof
- nonzero_tol -nonzero tolerance

## Outputs

- A -cleaned-up array

## Implementation structure

- Array clean-up utility. Drops non-zero elements with magnitude below the
- user-specified tolerance and converts between sparse and full storage de-
- pending on the density of non-zeroes in the array. Syntax:
- A=clean_up(spin_system,A,nonzero_tol)
- A -a numerical array or a cell array thereof
- nonzero_tol -nonzero tolerance
- A -cleaned-up array
- Skip opium objects
- Skip if disabled
- Process cells recursively
- Process polyadics recursively
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isnan()`, `iscell()`, `grumble()`, `ismember()`, `issparse()`, `nnz()`, `any()`, `all()`.
