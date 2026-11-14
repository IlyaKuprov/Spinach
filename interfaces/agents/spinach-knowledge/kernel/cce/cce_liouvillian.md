# kernel/cce/cce_liouvillian.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_liouvillian.m)

## Purpose

Build a local GKSL generator for full density or conditional coherence blocks.

## Behaviour

- Embeds explicitly supported jump operators in physical tensor order and uses column-major vectorization.
- Includes nonnegative rates in inverse seconds and local cluster Hamiltonians in rad/s; incomplete correlated jump support is refused.
- Disjoint exterior jumps are omitted only when they do not supply a frozen exterior field; evolving exterior mean fields are not replaced by a dissipator.

## Inputs and outputs

```matlab
[G,meta]=cce_liouvillian(model,sites,context,parameters)
```

**Inputs**

- `model`, `sites`, `context`, `parameters` — prepared cluster inputs and `dissipation.jumps` records with ordered support IDs, operator matrix and rate in s^-1.

**Outputs**

- `G` — dense GKSL generator acting on column-vectorized density/coherence blocks.
- `meta` — tensor order, dimensions, Hamiltonians and included jump records.

Conditional jumps must preserve both selected branches under the jump operator and its loss operator `J' * J`; otherwise use `megcce`, including for central registers larger than two levels.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_liouvillian.m)
- [cce_master.m](cce_master.md)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [Method reference](https://doi.org/10.1103/PhysRevLett.132.250401)
