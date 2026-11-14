# kernel/cce/cce_prepare.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_prepare.m)

## Purpose

Compile canonical Spinach metadata for local cluster Hamiltonians.

## Behaviour

- Applies the requested stock `assume` approximation without building a basis, and retains the source spin system for the reference backend.
- Builds stable physical-ID maps, multiplicities, coordinates in Angstrom, local linear/quadratic terms and ordered pair tensors in rad/s from metadata already converted by `create`.
- Uses `hamiltonian_backend=auto` unless explicitly set. The direct domain includes labframe, nmr, esr, deer and deer-zz subject to the source checks; unsupported direct requests are rejected by the Hamiltonian builder.

## Inputs and outputs

```matlab
model=cce_prepare(spin_system,parameters)
```

**Inputs**

- `spin_system` — output of `create`, not unconverted `sys`/`inter`.
- `parameters` — `central_spins`, `bath_spins`, `assumptions` (name or `{name,retention}`), Euler `orientation` in rad, `interaction_cutoff` in Angstrom, and `numerics.atol`, `rtol`, `max_memory_bytes`; optional `central` and `hamiltonian_backend`.

**Outputs**

- `model` — canonical source metadata, physical IDs and maps, local and pair interactions, central projectors when supplied, backend capability and content hash.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_prepare.m)
- [cce_hamiltonian.m](cce_hamiltonian.md)
- [cce_project.m](cce_project.md)
