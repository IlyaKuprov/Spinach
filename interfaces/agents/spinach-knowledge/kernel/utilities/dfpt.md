# kernel/utilities/dfpt.m

- Signature: `subgraphs=dfpt(conmatrix,max_sg_size)`

## Purpose

Analyzes a system connectivity graph by crawling in all available directions to produce connected subgraphs up to the user-specified size.

## Physical / mathematical content

- The connectivity matrix identifies connected spins; each output row identifies the spins belonging to a subgraph.

## Numerical / algorithmic content

- Starts with one subgraph per spin. At each growth step, it finds neighboring spins not already in each subgraph and extends the subgraph in every available direction. Subgraphs without neighbors retain a repeated index during growth. It then sorts the indices within each row and removes duplicate rows.
- Converts the resulting subgraph indices to a sparse logical matrix.

## Parameters / inputs

- conmatrix -[nspins x nspins] matrix with 1 for
- connected spins and 0 elsewhere
- max_sg_size -maximum connected subgraph size

## Outputs

- subgraphs -[n_subgraphs x nspins] matrix; each
- row contains 1 for spins that belong
- to the subgraph and 0 for spins that
- do not

## Implementation structure

- Checks that `conmatrix` is logical and square and that `max_sg_size` is a real positive integer.
- Performs growth steps for sizes `2` through `max_sg_size`, then constructs the sparse logical output.