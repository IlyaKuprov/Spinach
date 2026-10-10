# kernel/utilities/dfpt.m

## Purpose

Graph partitioning module. Analyzes the system connectivity graph and creates a list of all connected subgraphs of up to the user-specified size by crawling the graph in all available directions.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dfpt.m>

## Behaviour

- Syntax: `subgraphs=dfpt(conmatrix,max_sg_size)`.
- A consistency check (`grumble`) enforces that `conmatrix` is a logical square matrix and `max_sg_size` is a real positive integer scalar; violations raise errors.
- The crawl starts at each spin: `subgraphs` is initialised as `uint32(1:nspins)'`, where `nspins=size(conmatrix,2)`.
- For each subgraph size from 2 to `max_sg_size`, the function loops over the current subgraphs. For each subgraph, it finds spins reachable from it via `any(conmatrix(:,subgraph),2)`, excludes spins already in the subgraph, and grows the subgraph in every direction by appending each neighbour.
- Isolated subgraphs (no reachable neighbours) get a dummy index: `[subgraph subgraph(end)]`.
- Subgraphs with neighbours are grown as `[repmat(subgraph,[numel(neighbours) 1]) neighbours]`.
- After each growth pass, the grown set is merged with `cell2mat`, each row is sorted with `sort(subgraphs,2)`, and duplicates are removed with `unique(subgraphs,'rows')`.
- The final subgraph array is returned as a sparse logical matrix: a row index is built with `repmat` over `1:size(subgraphs,1)` repeated `max_sg_size` times and flattened, the subgraph entries are flattened with `subgraphs(:)`, and `sparse(row_index,subgraphs,1,nsg,nspins)` is converted with `logical`.

## Inputs and outputs

Inputs:

- `conmatrix` — `[nspins x nspins]` matrix with 1 for connected spins and 0 elsewhere; must be logical and square.
- `max_sg_size` — maximum connected subgraph size; must be a real positive integer scalar.

Outputs:

- `subgraphs` — `[n_subgraphs x nspins]` sparse logical matrix; each row contains 1 for spins that belong to the subgraph and 0 for spins that do not.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=dfpt.m>
