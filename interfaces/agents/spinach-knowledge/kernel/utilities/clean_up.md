# kernel/utilities/clean_up.m

- Signature: `A=clean_up(spin_system,A,nonzero_tol)`

## Purpose

Array clean-up utility. Drops non-zero elements with magnitude below the user-specified tolerance and converts between sparse and full storage de- pending on the density of non-zeroes in the array. Syntax: A=clean_up(spin_system,A,nonzero_tol)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
