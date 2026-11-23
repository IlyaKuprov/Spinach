# cce_interlace — native conditional state averaging

## Native role and invocation

`raw=cce_interlace(spin_system,cluster,clusters,ensemble,pulse_sequence,parameters)` is an internal service called by `cce_schedule` with the experiment function handle. Public calculations enter through `clexpn`; native interlaced controls are nested in `parameters.cce`. Every recursive conditional calculation is evaluated by `cce_evaluate`.

## Conditional averaging

Interlaced sampling changes selected bath assignments on an explicitly constructed conditioning domain while preserving the outer realisation elsewhere. The implementation forms maximal conditioning cylinders, removes contained cylinders, enumerates supported basis levels, and deduplicates assignment rows. Enumeration is bounded by `parameters.cce.interlaced.max_states` and the resolved memory budget; these limits are explicit rejection thresholds, not silent truncations.

Each assignment is passed to `cce_evaluate`. Conditional means retain absolute error and validity evidence. For an irreducible cluster, required subset responses are recomputed under the same conditioning domain before denominator formation. Error propagation, denominator-resolution checks, and invalid masks accompany the quotient; unresolved factors remain invalid. Per-state nonlinear reduction must precede the outer ensemble average.

With `raw_only=true`, the service returns the conditional physical raw average instead of a recursively divided factor. This is an engine control, not an additional public experiment. The returned `kind` distinguishes irreducible and exact full-bath records; `assignment_count` records conditional assignments.

## Native restrictions

`clexpn` rejects interlaced `noise_correlation`: connected noise uses ordinary ensemble moment averaging. Checkpoint ownership remains with the native scheduler; the interlaced service returns raw records rather than publishing checkpoint files itself. The experiment function handle, not a method-name dispatch, determines physical propagation.

[MATLAB source](../../../../../kernel/cce/cce_interlace.m)

## Native reference validation

Before enumerating a native conditional ensemble, the reference bath IDs and assignments must match the spin system. Every index must select a local basis vector and bath_rho must equal its projector within dimension^2*eps in Frobenius norm. Shape or density inconsistencies raise cce:spin_config before evaluation. This check does not use the propagation error tolerance; mixed reference densities belong to an explicit `bath_state`.

The conditioning contract uses cce.interlaced.max_states. Returned assignments, weights, condition_mass, log_condition_mass, and distribution_signature describe the deduplicated physical distribution; weights use the same shifted logarithmic product calculation as the established conditional average. raw.value, raw.valid, and raw.abs_error_bound always remain conditional raw quantities. factor, factor_valid, and factor_abs_error_bound separately carry the irreducible quotient; raw_only suppresses that quotient. kind identifies irreducible or raw_full_bath records. Every subset reuses the parent's cylinders. The allocation preflight includes assignment sorting workspace, and residual metadata retains the worst constituent solve, not merely the first assignment.

A direct native interlacing request rejects cce.checkpoint, for both writes and resumes, before constructing its conditioning distribution. Interlaced calculations do not own an ordinary raw-task checkpoint ledger.

All six arguments are mandatory; the service accepts native ensemble records and an experiment function handle only. Its local validator checks reference indices, required native records, and checkpoint ownership before conditional enumeration.

The local entry validator checks the state budget and the common reference density/assignment before recursive conditional evaluation. The recursion changes conditioning cylinders but retains that same validated reference.
