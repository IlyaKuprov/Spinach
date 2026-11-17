# kernel/cce/cce_liouvillian.m

## Native role and invocation

Expert local-jump Liouvillian construction service. Native retained generator construction is owned by cce_generators.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[G,meta]=cce_liouvillian(model,sites,spin_config,parameters)`. Engine records are not public experiment parameter structures.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_liouvillian.m)

## Purpose

Build a local GKSL generator for full density or projected signal blocks.

## Behaviour

Projected ME evolution checks invariance of both supplied observed manifold subspaces under the actual cluster Hamiltonian before constructing its generator. Nonzero observed-transition leakage raises `cce:population_transfer` and requests `megcce`; full generalized evolution is unchanged.

- The per-call memory ceiling includes full observed-spin system physical operators, the projected Liouvillian and both retained embedded matrices for every fully included jump, and is checked before constructing any of them.


- Projected manifold-leakage checks scale with the jump and loss-operator norms; reciprocal jump/rate rescaling cannot conceal manifold-changing dynamics.

- Embeds explicitly supported jump operators in physical tensor order and uses column-major vectorization.
- Includes nonnegative rates in inverse seconds and local cluster Hamiltonians in rad/s; incomplete correlated jump support is refused.
- Disjoint exterior jumps are omitted only when they do not supply a frozen exterior field; evolving exterior mean fields are not replaced by a dissipator.

## Inputs and outputs

```matlab
[G,meta]=cce_liouvillian(model,sites,spin_config,parameters)
```

**Inputs**

- `model`, `sites`, `spin_config`, `parameters` — prepared cluster inputs and `parameters.jumps` records with ordered support IDs, operator matrix and rate in s^-1.

**Outputs**

- `G` — dense GKSL generator acting on column-vectorized density/signal blocks.
- `meta` — tensor order, dimensions, Hamiltonians and included jump records.

Projected jumps must preserve both selected manifolds under the jump operator and its loss operator `J' * J`; otherwise use `megcce`, including for observed registers larger than two levels.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_liouvillian.m)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [Method reference](https://doi.org/10.1103/PhysRevLett.132.250401)

Validation is local to this service and its consumed numerical records. Public experiment and approximation options are validated by `cce_options` before native context dispatch.
