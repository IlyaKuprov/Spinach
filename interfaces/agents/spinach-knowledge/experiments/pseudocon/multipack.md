# experiments/pseudocon/multipack.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/multipack.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=multipack.m)

## Purpose

Splits a linear stream of spherical multipole moments into rank-specific cell entries. This is data packing: it does not fit a tensor, rotate moments, or run a pulse sequence.

## Inputs and output

- `ranks`: a real row vector of unique non-negative integers; each entry is a spherical rank L.
- `moments`: numeric data with exactly sum(2*L+1) elements across the declared ranks. Elements are consumed in linear order.
- `Ilm`: a cell array with the same dimensions as `ranks`. Each cell receives the next 2*L+1 values for that rank, in the input rank order.

The block sizes follow the multipole convention cited in the source: 1 component for L=0, 3 for L=1, 5 for L=2, and so on. The routine only partitions the stream; it imposes no ordering on ranks beyond preserving the supplied order.

## Checks and reference

The implementation rejects non-row, non-real, non-integer, negative, or repeated ranks, nonnumeric moments, and a moment count that differs from the sum of the rank block sizes. It does not separately require `moments` to be a vector; indexing is linear.

The source cites [DOI 10.1039/C6CP05437D](https://doi.org/10.1039/C6CP05437D) for the multipole-moment convention.