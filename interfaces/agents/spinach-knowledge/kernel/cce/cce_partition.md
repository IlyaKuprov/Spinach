# kernel/cce/cce_partition.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)

## Purpose

Freeze bounded physical spin groups for partitioned CCE.

## Behaviour

- Supports user groups, deterministic coupling/spatial grouping, and pilot-based adaptive regrouping; the resulting groups are frozen before production.
- Enforces physical size and Hilbert-dimension bounds without imposing graph-edge thresholds on strongly coupled members within a group.
- Adaptive proposals are accepted only if the declared pilot objective falls. Validation-only groups, pilot configuration/state IDs, and held-out production ranges prevent reusing pilot draws.

## Inputs and outputs

```matlab
partition=cce_partition(model,graph,parameters)
```

**Inputs**

- `model`, `graph` — canonical model and bath connectivity.
- `parameters.partition` — `mode`, `max_group_size`, `max_group_dimension`, and mode-specific `groups` or pilot fields. `order_semantics` is `max_groups` (default), `disjoint`, or `disjoint_inter`; see the source header for pilot score and holdout fields.

**Outputs**

- `partition` — frozen groups and IDs, group-level graph/order metadata, and pilot/validation diagnostics.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)
- [cce_clusters.m](cce_clusters.md)
- [cce_graph.m](cce_graph.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.110.L220302)
