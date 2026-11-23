# kernel/cce/cce_clusters.m

## Native role and invocation

Internal deterministic cluster-family enumerator. clexpn supplies the compiled spin system and the resolved parameters.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `clusters=cce_clusters(spin_system,graph,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_clusters.m)

## Purpose

Enumerate a deterministic subset-closed CCE cluster family.

## Behaviour

- When `sys.tols.cce_mem` is finite, graph/frontier growth, closure/sorting copies, decimal family signatures, and final cluster/subset metadata are checked before allocation (`cce_clusters:memory`). A cluster-count cap does not substitute for a byte budget; standalone combinatorial calls without `numerics` retain their existing unbudgeted contract.
- Enumerates connected candidate clusters through the requested order; disconnected subsets required by retained clusters are restored for recursive reduction.
- `periodic_segments` requires a single cycle and a truncated order and keeps only contiguous segments; disconnected subsets are not restored.
- Applies optional per-order cluster caps and explicit group-cluster selection deterministically, preserving physical-ID and subset bookkeeping.
- Partition mode expands group clusters into physical sites, so the reported order counts groups, not physical spins.

## Inputs and outputs

**Inputs**

- `spin_system`, `graph` — compiled spin system and physical connectivity.
- `parameters` — `order`, optional `groups` (default singleton sites), explicit group-ID `clusters`, and `cluster_cap` (default Inf); frozen partition groups are passed in by the driver.

**Outputs**

- `clusters` — ordered records with stable `id`, physical `sites`, `order`, `groups`, full `dimension`, `bath_dimension`, `proper_subsets`, `direct_children`, signed `mobius_power`, `connected_components`, `config_id`, `inclusion_reason`, selection `weight`, and exact `family_id`. `selected_counts`, `actual_counts`, and `cap_overrun` report the closure cost at each order.

## References

- [cce_graph.m](cce_graph.md)
- [cce_partition.m](cce_partition.md)
- [cce_reduce.m](cce_reduce.md)

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.

The sparse bath-to-group incidence matrix collapses physical edges and inter-group scores; within-group diagonals are excluded from ranking.
