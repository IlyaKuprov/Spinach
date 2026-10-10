# kernel/utilities/prune_subgraphs.m

## Purpose

Removes subgraphs that are contained entirely within other subgraphs, keeping only subgraphs that are not proper subsets of any other subgraph in the set.

## Behaviour

- The function first validates its input through an internal consistency check (`grumble`), which errors with `'subgraphs must be a logical array.'` if the input is not logical.
- Trivial cases are ignored and returned unchanged: when the input has fewer than 2 rows (`ngraphs < 2`) or fewer than 2 columns (`nspins < 2`), the function returns immediately.
- Spin counts per subgraph are computed as row sums: `spin_counts = sum(subgraphs,2)`.
- The subgraph overlap matrix is computed as `C = subgraphs*transpose(subgraphs)`, so `C(i,j)` is the number of spins shared between subgraphs `i` and `j`.
- A subgraph `i` is flagged as contained in subgraph `j` when `C(i,j)` equals the spin count of `i` (full containment of `i` in `j`) and subgraph `j` is strictly larger (`transpose(spin_counts) > spin_counts`, i.e. `spin_counts(j) > spin_counts(i)`).
- Pruning removes all rows flagged as subsets of any other subgraph: `subgraphs = subgraphs(~any(supersets,2),:)`.

## Inputs and outputs

**Syntax:**

```
subgraphs = prune_subgraphs(subgraphs)
```

**Input:**

- `subgraphs` — `[ngraphs x nspins]` logical array with `1` when a spin belongs to a subgraph and `0` otherwise. Must be logical; a non-logical input raises an error.

**Output:**

- `subgraphs` — `[ngraphs x nspins]` logical array with `1` when a spin belongs to a subgraph and `0` otherwise, with all subgraphs fully contained in a strictly larger subgraph removed.

## References

- Source: [prune_subgraphs.m on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/prune_subgraphs.m)
- Spinach Wiki: [prune_subgraphs.m](https://spindynamics.org/wiki/index.php?title=prune_subgraphs.m)
