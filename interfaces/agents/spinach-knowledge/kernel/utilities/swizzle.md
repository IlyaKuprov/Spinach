# kernel/utilities/swizzle.m

- Signature: `tuples=swizzle(index_arrays)`

## Purpose

Flattens out nested index lists and outputs them as an array of tuples in random order. This is useful for flattening nes- ted loops for parallel processing. Syntax: tuples=swizzle(index_arrays)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- index_arrays -a cell array of row vectors

## Outputs

- tuples -a matrix of tuples in random or-
- der, with tuples listed as rows

## Implementation structure

- Flattens out nested index lists and outputs them as an array
- of tuples in random order. This is useful for flattening nes-
- ted loops for parallel processing. Syntax:
- tuples=swizzle(index_arrays)
- index_arrays -a cell array of row vectors
- tuples -a matrix of tuples in random or-
- der, with tuples listed as rows
- Check consistency
- Kronecker up the arrays
- Randomise the tuple list
- Consistency enforcement
- Acording to a conference rumour, before appointing IK to a tenured
