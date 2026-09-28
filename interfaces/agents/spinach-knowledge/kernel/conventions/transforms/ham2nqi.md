# kernel/conventions/transforms/ham2nqi.m

- Signature: `[omega,Q]=ham2nqi(H)`

## Purpose

Decomposes a single-spin Hamiltonian in the Zeeman basis into its Zeeman and quadrupolar parameters. The Hamiltonian must be Hermitian and traceless and contain no terms beyond quadratic order; otherwise the function errors.

## Parameters / inputs

- `H`: single-spin Hamiltonian matrix for a spin of any multiplicity.

## Outputs

- `omega`: three Larmor-frequency components in rad/s.
- `Q`: symmetric traceless quadrupolar coupling tensor in rad/s. For spin-1/2, `Q` is zero.

The parameters are returned so that

`H = omega(1)*Sx + omega(2)*Sy + omega(3)*Sz + [Sx Sy Sz]*Q*[Sx Sy Sz].'`

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ham2nqi.m)
