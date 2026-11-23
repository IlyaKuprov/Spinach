# cce_compile — compiled cluster-expansion spin system

[MATLAB source](../../../../../kernel/cce/cce_compile.m)

## Native role and invocation

Internal compiler from create output and native parameters to the compiled spin system; no full-bath basis is constructed.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `spin_system=cce_compile(spin_system,parameters,assumptions)`. Engine records are not public experiment parameter structures.

Keywords: cce_compile clexpn canonical tensors physical IDs assumptions Zeeman formalism groups units

`cce_compile` accepts create output and parameters already resolved by `cce_options`; it does not resolve defaults itself. Tolerances are read from `spin_system.tols.cce_atol`, `cce_rtol` and `cce_mem`, which `create` fills. The formalism is `zeeman-liouv` when `parameters.bath_jumps` is nonempty; otherwise an existing `spin_system.bas.formalism` is kept, and `zeeman-hilb` is used when there is none. Bath IDs are parameters.cce.bath_spins, not cce.bath. An absent interaction cutoff retains the complete finite interaction model.

The returned spin system carries the requested assumptions, `bas.formalism`, and a `cce` subfield with canonical rad/s tensors (zeeman, quadratic, pairs), the resolved transition with its metadata and projectors, the assumption arguments, a roundoff bound, and a source hash. Observed and bath spins, backend, and orientation stay in `parameters`; multiplicities, isotopes, and coordinates are read from the stock `comp` and `inter` fields. It applies assume without building a whole-bath basis. Its source_hash and caller-owned subsystem cache signatures include assumptions and representation. It never converts canonical tensors a second time or resamples isotope occupancy. Occupied sites and their physical-ID groups are those already supplied by spin_bath/create.

An absent basis selects zeeman-hilb. An explicit zeeman-hilb or zeeman-liouv selects that complete local representation; global basis arrays are discarded, and other formalisms fail. The implementation ports the numerical compiler algebra directly and preserves native parameter ownership. Internal cce_hamiltonian can consume its retained tensor metadata. Context resource resolution and generator construction are separate services.
Observed transitions are resolved once after native tensor compilation. Coherent mixed-product models with product-ms selectors reuse the observed-only native generator; other representations and bath-state modes retain stock observed-Hamiltonian construction. No bath mean field is included in the isolated observed transition.

The physical source fingerprint excludes cce.checkpoint, cce.retain_raw, and cce.convergence controls. These fields remain in the resolved parameters for their respective owners; excluding them from the physical hash lets resumed and newly evaluated raw generator records retain identical physical provenance when only scheduling or resume switches change. Tensors, Spinach tolerances, time meshes, preparations, bath states, and all other physical options remain part of the fingerprint.

The compiler validates its consumed spin identifiers, orientation, assumptions, resolved numerical budgets, and representation locally before tensor compilation. cce_options supplies the shared native defaults and vocabulary; no compiler accuracy or memory policy is invented.

Compilation visits occupied unordered pairs in ascending physical-ID order. Each pair stores its unprojected laboratory-frame tensor for graph scoring and its independently masked tensor for propagation; cancelling bare terms need not cancel after masking.
