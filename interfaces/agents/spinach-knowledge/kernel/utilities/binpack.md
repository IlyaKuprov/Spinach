# kernel/utilities/binpack.m

- Signature: `bins=binpack(box_sizes,bin_size)`

## Purpose

A simple 1D bin packing algorithm. Collects the list of numbers supplied into sublists that sum to the number that is smaller or equal to the number specified. The algorithm is not optimal, but it does the job. Syntax: bins=binpack(box_sizes,bin_size)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- box_sizes -a row vector of box sizes
- bin_size -an integer specifying the bin size

## Outputs

- bins -a cell array of index vectors specifying
- boxes allocated into each bin

## Implementation structure

- A simple 1D bin packing algorithm. Collects the list of numbers
- supplied into sublists that sum to the number that is smaller or
- equal to the number specified. The algorithm is not optimal, but
- it does the job. Syntax:
- bins=binpack(box_sizes,bin_size)
- box_sizes -a row vector of box sizes
- bin_size -an integer specifying the bin size
- bins -a cell array of index vectors specifying
- boxes allocated into each bin
- Check consistency
- Number the boxes
- Find boxes that are bigger than bins
