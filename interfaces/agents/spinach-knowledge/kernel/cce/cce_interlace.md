# cce_interlace — native conditional state averaging

## Native role and invocation

`raw=cce_interlace(model,cluster,ensemble,schedule,parameters,evaluator)` is an internal service called by `cce_schedule` with a registered experiment evaluator. Public calculations enter through `clexpn`; native interlaced controls are nested in `parameters.cce`. The evaluator is propagated through every recursive conditional calculation.

## Conditional averaging

Interlaced sampling changes selected bath assignments on an explicitly constructed conditioning domain while preserving the outer realization elsewhere. The implementation forms maximal conditioning cylinders, removes contained cylinders, enumerates supported basis levels and deduplicates assignment rows. Enumeration is bounded by `max_combinations` and the resolved memory budget; these limits are explicit rejection thresholds, not silent truncations.

Each assignment is passed to the registered retained-cluster evaluator. Conditional means retain absolute error and validity evidence. For an irreducible cluster, required subset responses are recomputed under the same conditioning domain before denominator formation. Error propagation, denominator-resolution checks and invalid masks accompany the quotient; unresolved factors remain invalid. Per-state nonlinear reduction must precede the outer ensemble average.

With `raw_only=true`, the service returns the conditional physical raw average instead of a recursively divided factor. This is an engine control, not an additional public experiment. The returned metadata includes `interlaced`, `interlaced_irreducible` and the number of conditional assignments.

## Native restrictions

`clexpn` rejects interlaced `noise_correlation`: connected noise uses ordinary ensemble moment averaging. Checkpoint ownership remains with the native scheduler; the interlaced evaluator returns raw records rather than publishing checkpoint files itself. The registered evaluator, not a method-name dispatch or an inferred callback signature, determines physical propagation.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_interlace.m)

## Native reference validation

Before enumerating a native conditional ensemble, the reference bath IDs and assignments must match the model. Every index must select a local basis vector and bath_rho must equal its projector within dimension^2*eps in Frobenius norm. Shape or density inconsistencies raise cce:spin_config with the same message as the legacy route. This check does not use the propagation error tolerance; mixed reference densities belong to explicit_product.
