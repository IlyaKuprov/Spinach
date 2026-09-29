# kernel/utilities/scomponents.m

## Purpose

Computes the strongly connected components of a directed graph using David Gleich's implementation of Tarjan's algorithm.

## Behaviour

- Validates the input via an internal consistency check (`grumble`), which errors with `'the input must be a square logical matrix.'` if the input is not logical, not a matrix, or not square.
- Converts the adjacency matrix to compressed sparse row (CSR) form via `sparse2csr(sparse(A))`, returning row pointers `rp` and column indices `ci`.
- Runs an iterative (explicit-stack) form of Tarjan's algorithm over all nodes `sv = 1:n`, skipping nodes already assigned to a root (`root(v) > 0`).
- Maintains per-node arrays: `root` (current component root), `dt` (discovery times, incremented by a counter `t`), and `sci` (component labels, set to `-1` while a node is on the stack).
- Uses a call stack `rs` of size `2*n` storing (node, row-index) pairs, and a component stack `cs` of size `n`.
- When a node's root equals itself, all nodes on the component stack down to that node are popped and assigned the current component number `cn`, which is then incremented.
- Component numbering starts at `1` and increases in the order components are finalised.

## Inputs and outputs

**Input:**

- `A` — a logical square matrix with `1` (true) marking connected nodes in the graph.

**Output:**

- `sci` — a column vector of integers specifying the strongly connected component each graph node belongs to.

## References

- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/scomponents.m>
- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=scomponents.m>
- R. E. Tarjan, algorithm reference cited in the source: <http://dx.doi.org/10.1137/0201010>
