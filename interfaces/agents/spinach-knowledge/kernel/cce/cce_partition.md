# kernel/cce/cce_partition.m

## Native role and invocation

Internal partition construction and hybrid-group metadata service used by clexpn before cluster enumeration.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `partition=cce_partition(model,graph,parameters)`. The native call takes the experiment parameters with partition controls under `parameters.cce`; the context supplies the physical model and registered pulse-sequence handle.

[MATLAB source](../../../../../kernel/cce/cce_partition.m)

## Purpose

Freeze bounded physical spin groups for partitioned CCE.

## Behaviour

- `numerics.max_memory_bytes` bounds partition graph/candidate storage, resident pilot/held-out configurations and patch diagnostic vectors, and retained patch cluster records before allocation (`cce_partition:memory`). Pilot and held-out configurations remain caller-supplied arrays.
- Supports user groups, deterministic coupling/spatial grouping, and pilot-based adaptive regrouping; the resulting groups are frozen before production. Pilot and held-out patch evaluations omit production checkpoint options; a fresh or resumed production ledger remains owned by `cce_schedule`.
- Enforces physical size and Hilbert-dimension bounds without imposing graph-edge thresholds on strongly coupled members within a group.
- Adaptive proposals are accepted only if the declared pilot objective falls. Validation-only groups, pilot configuration/state IDs, and held-out production ranges prevent reusing pilot draws.

## Inputs and outputs

**Inputs**

- `model`, `graph` — canonical model and bath connectivity.
- `parameters.cce.partition` — `mode=user|coupling|spatial|adaptive` and `max_dimension` (group Hilbert-dimension ceiling). Spatial grouping also requires `group_count` and `max_rounds`. Adaptive mode or optional `adaptive=true` requires `patches`, `patch_atol`, `patch_rtol`, `pilot_configs`, `heldout_configs`, `max_rounds`, and `max_probes`.
- `parameters.cce.groups`, `parameters.cce.expansion`, and `parameters.cce.order` — explicit groups and expansion-order semantics are CCE options, not partition fields.

**Outputs**

- `partition` — frozen groups and IDs, group-level graph/order metadata, and pilot/validation diagnostics.

Adaptive proper-bath patches restrict explicit adjacency, supplied scores, and resonance mismatch matrices to the patch physical-ID order before building the local graph.

Adaptive pilot and held-out validation use the raw complete-patch record when that patch is present in the selected family (or the raw baseline for an empty patch). Singular auxiliary ratios therefore cannot invalidate an otherwise valid exact patch, including generalised density matrices with zero entries and either supported zero policy.

## References

- [cce_clusters.m](cce_clusters.md)
- [cce_graph.m](cce_graph.md)
- [Method reference](https://doi.org/10.1103/PhysRevB.110.L220302)

Unresolved (infinite-objective) reductions may take a deterministic capacity-limited pair merge even when the objective remains infinite. Each merge reduces the number of groups, and round/probe ceilings still apply. Such an intermediate is never accepted: the selected partition must pass both pilot and held-out patch checks, including the exact-patch bypass when a full patch becomes reachable.

## Native adaptive partition contract

The context passes the native parameters structure directly to cce_partition(model,graph,parameters). Controls live in parameters.cce.partition; pilot_configs and heldout_configs contain native bath-state records (bath_ids, bath_rho, assignment, configuration_id, state_id, and weight). The context supplies model.pulse_sequence. Pilot and held-out patch calculations use that registered experiment evaluator and native reduction, preserving the declared time grid, finite controls, exterior state, and observable. Partition selection uses a bounded merge search, memory preflight, independent held-out validation, and completion diagnostics.

Native adaptive patch selection removes the experiment-owned propagator preference only from its internal graph and cluster-selection request. The native evaluator still receives the original numerical preference, so an explicit stock master-equation propagator is used unchanged for pilot, held-out, and production dynamics.

Only the native experiment parameter contract is accepted. Exact-patch evaluation always obtains its callback from the registered pulse sequence; no method-name dispatch or event-schedule compiler is used.

Native partition bounds and group membership are checked before constructing the internal graph and cluster-selection request.

Coupling merges consider only groups joined by a nonzero physical edge. Spatial swaps and adaptive pilot proposals may involve uncoupled groups.
