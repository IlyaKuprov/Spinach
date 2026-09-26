# kernel/utilities/iseye.m

- Signature: `verdict=iseye(M)`

## Purpose

Returns true for unit matrices. The test is designed to be computationally affordable. Syntax: verdict=iseye(M)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- M -a matrix

## Outputs

- verdict -true or false

## Implementation structure

- Returns true for unit matrices. The test is designed to be
- computationally affordable. Syntax:
- verdict=iseye(M)
- M -a matrix
- verdict -true or false
- Check consistency
- Run the checks
- Not even square
- Not even diagonal
- Test vector
- Compare with unit
- Test failed
