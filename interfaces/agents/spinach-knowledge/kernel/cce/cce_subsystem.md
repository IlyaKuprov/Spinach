# cce_subsystem — retained physical-spin ownership

[MATLAB source](../../../../../kernel/cce/cce_subsystem.m)

## Native role and invocation

Internal physical-ID-preserving stock subsystem and basis constructor used only after bath-cluster selection.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[spin_system,cache]=cce_subsystem(spin_system,parameters,sites,configuration,cache)`. Engine records are not public experiment parameter structures.

Keywords: cce_subsystem clexpn physical IDs ordering cache density meanfield isotopes multiplicities

`cce_subsystem` accepts cce_compile output, the resolved parameters, an ordered row of distinct physical bath IDs, a scalar configuration, and [] or a previously returned cache. Empty sites retains the full observed subsystem. Tensor-product order is [parameters.observed sites], never sorted physical IDs.

The returned spin system carries, in its `cce` subfield, sites, ids, physical_to_local (zero for discarded spins), Hilbert/observed/bath dimensions, retained groups, local Pauli matrices, bath_rho, and a bound configuration. It contains no stale whole-system basis arrays and allocates no global bath state or operator embeddings.

Explicit configurations give bath_ids and bath_rho for the entire selected bath, in matching order. Optional expectations must agree with those densities. Omitted densities use `parameters.bath_state`, or the maximally mixed state when it is absent; sampled and interlaced expansions require a supplied realisation. Local densities are positive, Hermitian, and trace one. The static exterior field is derived from the same density assignment for every cluster, and explicit exterior IDs/expectations must agree with it. With meanfield='none' the bound exterior arrays are empty. configuration_id, state_id, weight, and backend are optional configuration metadata.

The explicit single-entry cache stores an exact spin-system/parameters/ordered-sites/configuration signature, the returned spin system, and hit/miss counters. Changes to tensors, assumptions, formalism, site order, or configuration invalidate it. There is no persistent cache. Generator services own their later numerical caches.

For static mean fields, a finite parameters.cce.meanfield_radius requires finite coordinates for observed and exterior spins once the product configuration is bound. Retain exterior spins whose minimum distance to any observed spin is at most the radius in Angstrom. Included-cluster geometry is irrelevant to this selection. An absent or infinite radius remains unbounded, and compilation alone does not require geometry.

Explicit configuration density shapes, positivity, and normalisation are checked at entry before consulting the cache. Agreement of supplied expectations and exterior fields with computed density expectations remains a computed-domain consistency test.
