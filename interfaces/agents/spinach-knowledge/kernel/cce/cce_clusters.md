# kernel/cce/cce_clusters.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_clusters.m)

## Purpose

Enumerate a deterministic subset-closed CCE cluster family.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Enumerates connected candidate clusters through the requested order; disconnected subsets required by retained clusters are restored for recursive reduction.
- Applies optional per-order cluster caps and explicit group-cluster selection deterministically, preserving physical-ID and subset bookkeeping.
- Partition mode expands group clusters into physical sites, so the reported order counts groups, not physical spins.

## Inputs and outputs

```matlab
clusters=cce_clusters(model,graph,parameters)
```

**Inputs**

- `model`, `graph` — prepared model and physical connectivity.
- `parameters` — `order`, optional `groups` (default singleton sites), explicit group-ID `clusters`, and `cluster_cap` (default Inf); frozen partition groups are passed in by the driver.

**Outputs**

- `clusters` — ordered records with stable IDs, physical `sites`, `order`, `groups`, `hilbert_dim`, proper-subcluster references, per-order reduction powers, score, tuple_id and closure_only.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_clusters.m)
- [cce_graph.m](cce_graph.md)
- [cce_partition.m](cce_partition.md)
- [cce_reduce.m](cce_reduce.md)
