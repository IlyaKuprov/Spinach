# kernel/utilities/sorensen.m

- Signature: `b=sorensen(rho_init,rho_targ)`

## Purpose

Computes the Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators (Equation 186 of https://doi.org/10.1016/0079-6565(89)80006-8).

## Physical / mathematical content

The function diagonalises the initial and target matrices, sorts their eigenvalues in the same order, and computes

`b = (sigma_init' * sigma_targ) / trace(rho_targ^2)`.

This is an exact unitary bound. The amount reachable with realistically available instrumental controls may be smaller; see https://doi.org/10.1080/00268979909483117.

## Parameters / inputs

- `rho_init` — initial density matrix, Hilbert space.
- `rho_targ` — target density matrix, Hilbert space.

## Output

- `b` — Sorensen bound.

## Implementation structure

The inputs must be numeric Hermitian matrices of the same size. The function checks these conditions, diagonalises both matrices, sorts their eigenvalues, and evaluates the bound above.

Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators. Equation 186 from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax: b=sorensen(rho_init,rho_targ)

- from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax:

ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=sorensen.m>