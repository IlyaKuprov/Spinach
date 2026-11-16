# kernel/cce/cce_partition.m

## Native role and invocation

Internal partition construction and hybrid-group metadata service used by cce_clusters.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `partition=cce_partition(model,graph,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)

## Purpose

Freeze bounded physical spin groups for partitioned CCE.

## Behaviour

- `numerics.max_memory_bytes` bounds partition graph/candidate storage, resident pilot/held-out configurations and patch diagnostic vectors, and retained patch cluster records before allocation (`cce_partition:memory`). Pilot and held-out configurations remain caller-supplied arrays.


- Supports user groups, deterministic coupling/spatial grouping, and pilot-based adaptive regrouping; the resulting groups are frozen before production. Pilot and held-out patch evaluations omit production checkpoint options; a fresh or resumed production ledger remains owned by `cce_schedule`.
- Enforces physical size and Hilbert-dimension bounds without imposing graph-edge thresholds on strongly coupled members within a group.
- Adaptive proposals are accepted only if the declared pilot objective falls. Validation-only groups, pilot configuration/state IDs, and held-out production ranges prevent reusing pilot draws.

## Inputs and outputs

```matlab
partition=cce_partition(model,graph,parameters)
```

**Inputs**

- `model`, `graph` — canonical model and bath connectivity.
- `parameters.partition` — `mode=user|coupling|spatial|adaptive` and `max_dimension` (group Hilbert-dimension ceiling). Spatial grouping also requires `group_count` and `max_rounds`. Adaptive mode or optional `adaptive=true` requires `patches`, `patch_atol`, `patch_rtol`, `pilot_configs`, `heldout_configs`, `max_rounds`, and `max_probes`.
- `parameters.groups`, `parameters.expansion`, and `parameters.order` — explicit groups and expansion-order semantics are top-level parameters, not partition fields.

**Outputs**

- `partition` — frozen groups and IDs, group-level graph/order metadata, and pilot/validation diagnostics.

Adaptive proper-bath patches restrict explicit adjacency, supplied scores, and resonance mismatch matrices to the patch physical-ID order before building the local graph.

Adaptive pilot and held-out validation use the raw complete-patch record when that patch is present in the selected family (or the raw baseline for an empty patch). Singular auxiliary ratios therefore cannot invalidate an otherwise valid exact patch, including generalized density matrices with zero entries and either supported zero policy.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_partition.m)
- [cce_clusters.m](cce_clusters.md)
- [cce_graph.m](cce_graph.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.110.L220302)

Unresolved (infinite-objective) reductions may take a deterministic capacity-limited pair merge even when the objective remains infinite. Each merge reduces the number of groups, and round/probe ceilings still apply. Such an intermediate is never accepted: the selected partition must pass both pilot and held-out patch checks, including the exact-patch bypass when a full patch becomes reachable.
