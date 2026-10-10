# kernel/utilities/sphten2zeeman.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sphten2zeeman.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sphten2zeeman.m)

## Purpose

Returns a projector matrix `P` that converts state vectors written in the spherical tensor basis set used by Spinach into state vectors written in the Zeeman basis set in Liouville space, via `rho_zeeman = P * rho_sphten`.

## Behaviour

- The function first calls an internal consistency check (`grumble`) that errors with `'this function is only available for sphten-liouv formalism.'` unless `spin_system.bas.formalism` is `'sphten-liouv'`.
- Each substance is converted independently from `bas.basis{n}`, using the multiplicities of `chem.parts{n}`. The resulting sparse matrices are assembled as a direct sum: no inter-substance coherences are introduced.
- Each local tensor product is divided by its Frobenius norm and multiplied by `sqrt(D_n)`, where `D_n` is the local Hilbert dimension. Consequently the unit coordinate maps to `vec(I_D_n)`, and its value equals the Hilbert trace divided by `D_n`.
- To convert concentration-weighted spherical-tensor states into physical trace-equals-concentration Zeeman vectors, divide each destination block of P by its local Hilbert dimension D_n. The returned P retains the stock operator-normalisation convention.
- The projector need not be square and may be huge.

## Inputs and outputs

**Inputs**

- `spin_system` — main Spinach data structure using the `sphten-liouv` formalism and including basis set information.

**Outputs**

- `P` — projector matrix used as `rho_zeeman = P * rho_sphten`.

## References

- Spinach Wiki: [sphten2zeeman.m](https://spindynamics.org/wiki/index.php?title=sphten2zeeman.m)
