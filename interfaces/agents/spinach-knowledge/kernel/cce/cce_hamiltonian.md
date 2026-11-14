# kernel/cce/cce_hamiltonian.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_hamiltonian.m)

## Purpose

Build one cluster Hamiltonian in central-register-then-bath tensor order.

## Behaviour

- Includes every central spin and the requested physical bath sites, including the empty-bath central baseline.
- Uses direct Cartesian operator assembly when supported; the stock path uses `kill_spin`, `basis`, `assume`, `hamiltonian` and `orientation` as a reference/fallback.
- Applies exterior product-state mean fields from the supplied context and enforces the declared memory limit. A context backend override supports direct/reference comparisons.

## Inputs and outputs

```matlab
[H,meta]=cce_hamiltonian(model,sites,context)
```

**Inputs**

- `model` — output of `cce_prepare`.
- `sites` — physical bath IDs; `[]` selects the central baseline.
- `context` — scalar structure with optional backend override, exterior_ids and exterior_expectations (one real Ix/Iy/Iz row per exterior ID). Exterior fields use retained model pair interactions.

**Outputs**

- `H` — Hermitian cluster Hamiltonian in rad/s.
- `meta` — physical tensor order, dimensions, local operators and backend metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_hamiltonian.m)
- [cce_prepare.m](cce_prepare.md)
- [cce_states.m](cce_states.md)
