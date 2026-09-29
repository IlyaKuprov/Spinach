# kernel/utilities/sphten2zeeman.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sphten2zeeman.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sphten2zeeman.m)

## Purpose

Returns a projector matrix `P` that converts state vectors written in the spherical tensor basis set used by Spinach into state vectors written in the Zeeman basis set in Liouville space, via `rho_zeeman = P * rho_sphten`.

## Behaviour

- The function first calls an internal consistency check (`grumble`) that errors with `'this function is only available for sphten-liouv formalism.'` unless `spin_system.bas.formalism` is `'sphten-liouv'`.
- The projector is preallocated as a sparse matrix with `prod(spin_system.comp.mults.^2)` rows and `size(spin_system.bas.basis,1)` columns, initially with zero nonzeros, using `spalloc`.
- The destination (Zeeman) basis is not normalised; a destination normalisation factor `destin_norm = sqrt(prod(spin_system.comp.mults))` is computed once.
- A `parfor` loop runs over the rows of `spin_system.bas.basis` (the source basis set). For each basis element:
  - The state `rho` is built as a Kronecker product over all spins `k`, using the irreducible spherical tensors `irr_sph_ten(spin_system.comp.mults(k))` selected by the basis-set index `spin_system.bas.basis(n,k)+1`.
  - The source basis is not normalised; a per-column source normalisation `source_norm = norm(rho(:),2)` is computed.
  - The column of the projector is written as `P(:,n) = destin_norm * rho(:) / source_norm`.
- The projector need not be square and may be huge.

## Inputs and outputs

**Inputs**

- `spin_system` — main Spinach data structure using the `sphten-liouv` formalism and including basis set information.

**Outputs**

- `P` — projector matrix used as `rho_zeeman = P * rho_sphten`.

## References

- Spinach Wiki: [sphten2zeeman.m](https://spindynamics.org/wiki/index.php?title=sphten2zeeman.m)
