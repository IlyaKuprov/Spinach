# kernel/utilities/sphten2zeeman.m

- Signature: `P=sphten2zeeman(spin_system)`

## Purpose

Constructs a projector that converts state vectors from Spinach’s spherical-tensor basis to the Zeeman basis in Liouville space. For a state vector `rho_sphten`, the corresponding vector is `rho_zeeman = P * rho_sphten`.

## Mathematical content

- For each source basis row, the function forms the Kronecker product of the corresponding per-spin irreducible spherical-tensor matrices, divides the resulting vectorized matrix by its 2-norm, and scales it by `sqrt(prod(spin_system.comp.mults))`. These factors account for the source and destination basis normalizations.
- `P` is sparse, with `prod(spin_system.comp.mults.^2)` rows and one column per row of `spin_system.bas.basis`. It need not be square.
- The input must use the `sphten-liouv` formalism.

## Parameters / inputs

- `spin_system` — Spinach data structure using the `sphten-liouv` formalism and containing basis-set information.

## Outputs

- `P` — projector matrix mapping spherical-tensor-basis state vectors to Zeeman-basis state vectors: `rho_zeeman = P * rho_sphten`. It may be large and need not be square.

## Implementation structure

- Preallocates a sparse projector with one column for each basis-set row.
- Builds each column from the tensor product of per-spin spherical-tensor matrices, applies the normalization factors, and fills the columns in a `parfor` loop.
- Rejects spin systems whose formalism is not `sphten-liouv`.

## References

- [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=sphten2zeeman.m)
