# kernel/cce/cce_partition.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)

## Purpose

Freeze bounded physical spin groups for partitioned CCE.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Supports user groups, deterministic coupling/spatial grouping, and pilot-based adaptive regrouping; the resulting groups are frozen before production.
- Enforces physical size and Hilbert-dimension bounds without imposing graph-edge thresholds on strongly coupled members within a group.
- Adaptive proposals are accepted only if the declared pilot objective falls. Validation-only groups, pilot configuration/state IDs, and held-out production ranges prevent reusing pilot draws.

## Inputs and outputs

```matlab
partition=cce_partition(model,graph,parameters)
```

**Inputs**

- `model`, `graph` — canonical model and bath connectivity.
- `parameters.partition` — `mode=user|coupling|spatial|adaptive` and `max_dimension` (group Hilbert-dimension ceiling). Spatial grouping also requires `group_count` and `max_rounds`. Adaptive mode or optional `adaptive=true` requires `patches`, `patch_atol`, `patch_rtol`, `pilot_contexts`, `heldout_contexts`, `max_rounds`, and `max_probes`.
- `parameters.groups`, `parameters.expansion`, and `parameters.order` — explicit groups and expansion-order semantics are top-level parameters, not partition fields.

**Outputs**

- `partition` — frozen groups and IDs, group-level graph/order metadata, and pilot/validation diagnostics.

Adaptive proper-bath patches restrict explicit adjacency, supplied scores, and resonance mismatch matrices to the patch physical-ID order before building the local graph.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)
- [cce_clusters.m](cce_clusters.md)
- [cce_graph.m](cce_graph.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.110.L220302)
