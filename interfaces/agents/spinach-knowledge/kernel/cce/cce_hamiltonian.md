# kernel/cce/cce_hamiltonian.m

## Native role and invocation

Expert retained-Hilbert Hamiltonian assembly service. The native context generator path is cce_subsystem followed by cce_generators, which can call this helper for validated coherent Hilbert assembly; it is not an alternative public context.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[H,meta]=cce_hamiltonian(spin_system,parameters,sites,spin_config)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_hamiltonian.m)

## Purpose

Build one cluster Hamiltonian in observed-spin system-then-bath tensor order.

## Behaviour

- Includes every observed spin and the requested physical bath sites, including the empty-bath observed baseline.
- Uses direct Cartesian operator assembly when supported; the stock path uses `kill_spin`, `basis`, `assume`, `hamiltonian`, and `orientation` as a reference/fallback.
- Applies exterior product-state mean fields from the supplied spin_config and enforces the declared memory limit. A spin_config backend override supports direct/reference comparisons.

## Inputs and outputs

**Inputs**

- `spin_system` — output of `cce_compile`.
- `parameters` — parameters resolved by `cce_options`.
- `sites` — physical bath IDs; `[]` selects the observed baseline.
- `spin_config` — scalar structure with optional backend override, exterior_ids, and exterior_expectations (one real Ix/Iy/Iz row per exterior ID). Exterior fields use retained pair interactions.

**Outputs**

- `H` — Hermitian cluster Hamiltonian in rad/s.
- `meta` — physical tensor order, dimensions, local operators, and backend metadata.

## References

- [cce_compile.m](cce_compile.md)

The retained service uses the native compilation/input boundary: physical models come from cce_compile, and generated baths enter canonical sys/inter through spin_bath.

Each unordered pair supplies its retained tensor in ascending site order and its transpose in reverse order. Missing pairs contribute no field.
